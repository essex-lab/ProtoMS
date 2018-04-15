# Author: Chris Cave-Ayland

from collections import defaultdict
from functools import partial
import glob
import matplotlib
import matplotlib.pyplot as pl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from operator import add
import os
import pymbar

import simulationobjects as sim
import free_energy_base as feb
import calc_ti_decomposed
import calc_gci2


if "DISPLAY" not in os.environ or os.environ["DISPLAY"] == "":
    matplotlib.use('Agg')

# logger = logging.getLogger('protoms')
thermal_wavelength = 1.00778365325  # of water, in angstroms
kB = 0.0019872  # Boltzmann's constant, in kcal.mol^-1.K^-1
standard_volume = 30.  # Standard volume of a water in bulk water
tip4p_excess = -6.2  # Excess chemical potential of tip4p in bulk tip4p


class PMF2D(feb.Series):
    def __init__(self, lambdas, Bs, volume, temperature, *args):
        self.coordinate = lambdas
        self.coordinate2 = Bs
        self.volume = volume
        self.temperature = temperature
        self.values = []
        for dat in zip(*args):
            self.values.append(feb.PMF(self.coordinate2, *dat))

    def equilibrium_B(self):
        """Get the value of B equating to equilibrium with bulk water."""
        beta = 1 / (self.temperature * kB)
        return beta * tip4p_excess + np.log(self.volume / standard_volume)

    @property
    def dG(self):
        eqb_B = self.equilibrium_B()
        closest = np.argmin([abs(B - eqb_B) for B in self.coordinate2])
        return self.values[-1].values[closest] - self.values[0].values[closest]

    def plot(self, axes, show_error=True,  **kwargs):
        lambdas, Bs = np.meshgrid(self.coordinate2, self.coordinate)
        values = np.array([[col.value for col in row] for row in self.values])
        axes.plot_surface(lambdas, Bs, values, cstride=1, rstride=1)
        return axes


class Result2D(feb.Result):
    @property
    def pmf(self):
        pmfs = [PMF2D(self.lambdas, dat[0].coordinate2,
                      dat[0].volume, dat[0].temperature, *dat)
                for dat in self.data]
        return reduce(add, pmfs)


class GCMCMBAR(feb.MBAR):
    results_name = 'results_inst'
    # subdir_glob = 'b_*/'
    pmf_class = PMF2D
    result_class = Result2D

    def __init__(self, lambdas, volume=30., **kwargs):
        self.lambdas = lambdas
        self.volume = volume
        self.data = []
        self._data_lambdas = []
        self.subdir_glob = 'b_*/'
        self._data_Bs = []

    def add_data(self, series):
        dat = np.array([series.feenergies[lam]
                        for lam in sorted(series.feenergies)])
        N = dat.shape[-1]
        # add solventson data as the final row of each data entry
        # not an ideal solution but saves the need to rewrite __getitem__
        dat = np.append(dat, series.solventson.reshape((1, N)), axis=0)
        self.data.append(dat)
        self._data_lambdas.append(series.lam[0])
        self._data_Bs.append(series.bvalue)
        return N

    def B_to_chemical_potential(self, B, temperature):
        """Convert a B value to the equivalent chemical potential"""
        beta = 1 / (temperature * kB)
        return (B-np.log(self.volume/thermal_wavelength**3)) / beta

    @property
    def Bs(self):
        return sorted(set(self._data_Bs))

    def calculate(self, temp=300.):
        N = self.data[0].shape[-1]
        N_sims = len(self.data)
        N_Bs = len(self.Bs)
        N_lams = len(self.lambdas)
        beta = 1 / (sim.boltz * temp)

        u_kn = np.zeros(shape=(N_sims, N*N_sims))
        ns = np.zeros(N*N_sims)
        for dat, B_sim, lam_sim in zip(
                self.data, self._data_Bs, self._data_lambdas):
            # due to globbing simulation data will not have been
            # read in a logical order, need the below indices to
            # determine the proper position for each set of data
            i_B = self.Bs.index(B_sim)
            i_lam = self.lambdas.index(lam_sim)
            ns = dat[-1]  # water occupancy  data
            es = dat[:-1]  # simulation energies
            for i in xrange(N_lams):
                for j, B in enumerate(self.Bs):
                    mu = self.B_to_chemical_potential(B, temp)
                    u_kn[i*N_Bs+j, (i_lam*N_Bs+i_B)*N: (i_lam*N_Bs+i_B+1)*N] =\
                        beta*(es[i] + ns * mu)
        samples = np.array([N] * N_sims)
        mbar = pymbar.MBAR(u_kn, samples)
        free_energies = mbar.getFreeEnergyDifferences(
            compute_uncertainty=False, return_theta=False)[0] / beta

        # free energy order seems to be reversed with respect to B values
        return PMF2D(self.lambdas, self.Bs[::-1], self.volume, temp,
                     free_energies[0].reshape([N_lams, N_Bs]))


class GCMC_MBAR(object):
    """Calculate free energy differences from simulations over
    multiple lambda and B values.

    Parameters
    ----------
    path : string
      folder_location
    file_name : string
      name of results files to load, should usually be results_inst
    temperature : float, kelvin
      simulation temperature
    volume : float, angstroms^3
      volume of gcmc simulation region
    skip : int
      number of configurations to skip from each results file
    max_read : int
      maximum number of configurations to read from each results file
    """
    def __init__(self, path, file_name, temperature, volume, skip, max_read):
        self.path = path
        self.file_name = file_name
        self.skip = skip
        self.beta = 1. / (temperature * kB)  # thermodynamic beta
        self.volume = volume
        self.max_read = max_read

    def _parse_folder(self, path, lambdas, Bs):
        """Parse ProtoMS results file in path and calculate the
        reduced potential energies for all provided lambda and B values.

        Parameters
        ----------
        path : string
          folder location
        lambdas : list of floats
          all lambda values to include in MBAR calculation
        Bs : list of floats
          all Bs to include in MBAR calculation

        Returns
        -------
        state_reduced_energies : numpy array of floats
          reduced potential energy evaluated at all states
        """

        results_file = sim.ResultsFile()
        results_file.read(filename=os.path.join(path, self.file_name),
                          skip=self.skip, readmax=self.max_read)
        results_file.make_series()
        state_reduced_energies = np.zeros((len(lambdas) * len(Bs),
                                           len(results_file.snapshots)))

        i = 0
        for lam in lambdas:
            for B in Bs:
                chem_pot = self.B_to_chemical_potential(B)
                energies = results_file.series.feenergies[lam]
                Es = self.beta * (energies + chem_pot *
                                  results_file.series.solventson)
                state_reduced_energies[i] = Es
                i += 1

        return state_reduced_energies

    def _calc_reduced_energy_matrix(self):
        """Get the reduced potential energy matrix required for MBAR

        Returns
        -------
        lambdas : list of floats
          sorted lambda values from results
        Bs : list of floats
          sorted B values from results
        reduced_energies : 3d numpy array of floats
          reduced energies for each configuration of each state
          evaluated for all states, suitable for use with MBAR
        """

        # Find lambdas and sort them
        paths = glob.glob(os.path.join(self.path, "lam-*"))
        paths.sort()
        gcfldr = glob.glob(paths[0] + "/b_*")

        lambdas = [float(p.split('lam-')[-1]) for p in paths]
        Bs = [float(p.split('b_')[-1]) for p in gcfldr]
        lambdas.sort()
        Bs.sort()

        reduced_energies = []
        for lamfldr in paths:
            gcmcfolders = glob.glob(lamfldr + "/b_*")
            gcmcfolders.sort(key=lambda x: float(x.split('/b_')[1]))
            for gcfldr in gcmcfolders:
                reduced_energies.append(self._parse_folder(
                    gcfldr, lambdas, Bs))

        reduced_energies = np.array(reduced_energies)
        reduced_energies.dump('old.pkl')
        return lambdas, Bs, reduced_energies

    def B_to_chemical_potential(self, B):
        """Convert a B value to the equivalent chemical potential"""
        return (B-np.log(self.volume/thermal_wavelength**3)) / self.beta

    def equilibrium_B(self):
        """Get the value of B equating to equilibrium with bulk water."""
        return self.beta * tip4p_excess + np.log(self.volume / standard_volume)

    def get_free_energies(self):
        """
        Evaluate free energy differences for simulations of multiple
        lambda and B values using MBAR

        Returns
        -------
        numpy array
          all lambda values
        numpy array
          all B values
        numpy array
          free energy differences for all states

        """

        # get the reduced energies
        lambdas, Bs, reduced_energies = self._calc_reduced_energy_matrix()

        samples = np.array([reduced_energies.shape[2]] *
                           reduced_energies.shape[0])
        mbar = pymbar.MBAR(reduced_energies, samples)
        free_energies = mbar.getFreeEnergyDifferences(
            compute_uncertainty=False, return_theta=False)[0] / self.beta
        free_energies = free_energies[0].reshape((len(lambdas), len(Bs)))

        return lambdas, Bs, free_energies


# def get_arg_parser():
#     import argparse
#     # Setup a parser of the command-line arguments
#     parser = argparse.ArgumentParser(
#         description="Program to calculate free energy from thermodynamic integration")
#     parser.add_argument('-d', '--directory', help="the root directory that contains all the output files of the simulation. Default is cwd.", default="./")
#     parser.add_argument('-r', '--results', help="the name of the file to analyse. Default is results. ", default="results_inst")
#     parser.add_argument('-s', '--skip', type=int, help="the number of blocks to skip to calculate the free energy differences in one window. default is 0. Skip must be greater or equal to 0", default=0)
#     parser.add_argument('-m', '--max', type=int, help="the upper block to use. default is 99999 which should make sure you will use all the available blocks. max must be greater or equal to 0", default=99999)
#     parser.add_argument('--plot', action='store_true', default=True, help="3D Plots the gradient and free energy surface against lambda and B values")
#     parser.add_argument('-t', '--temperature', default=298.15, help="Simulation temperature.")
#     parser.add_argument('-v','--volume', type=float, help="volume of the GCMC insertion region", required=True)
#     return parser

def plot_surface(lambdas, Bs, Z):
    fig = pl.figure()
    ax = Axes3D(fig)
    lambdas, Bs = np.meshgrid(Bs, lambdas)

    ax.plot_surface(lambdas, Bs, Z,
                    cstride=1, rstride=1)
    return fig


# if __name__ == '__main__':
#     args = get_arg_parser().parse_args()

#     # Fix negative values of skip and max
#     if args.max < 0:
#         args.max = 99999
#     if args.skip <= 0:
#         args.skip = -1

#     mbar = GCMC_MBAR(args.directory, args.results, args.temperature,
#                      args.volume, args.skip, args.max)
#     lambdas, Bs, free_energies = mbar.get_free_energies()

#     eqb_B = mbar.equilibrium_B()
#     closest = np.argmin([abs(B - eqb_B) for B in Bs])

#     print "Equilibrium B value is %.4f." % eqb_B
#     print "Closest simulated B value is %.4f." % Bs[closest]
#     print "Alchemical free energy difference at B=%.4f is %.4f" % (
#         Bs[closest], free_energies[-1][closest] - free_energies[0][closest])

#     if args.plot:
#         fig = plot_surface(lambdas, Bs, free_energies)
#         fig.savefig('free_energy_surface.png')
#         pl.show()
