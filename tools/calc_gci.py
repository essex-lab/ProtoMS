# Author: Gregory Ross
"""
Routines to calculate free energies with GCMC as outlined in
J. Am. Chem. Soc., 2015, 137 (47), pp 14930-14943

Can be executed from the command line as a stand-alone program
"""

import sys
import numpy as np
import os

import free_energy_base as feb
from gcmc_free_energy_base import GCI, GCMCResult, tip4p_excess
from table import Table

import matplotlib
if "DISPLAY" not in os.environ or os.environ["DISPLAY"] == "":
    matplotlib.use('Agg')
import matplotlib.pyplot as plt


class TitrationCalculation(feb.FreeEnergyCalculation):
    """Calculate free energy differences using the Grand Canonical Integration
    method. Collates water occupancy data from simulations carried out at
    diferent chemical potentials and gives NVT binding free energies."""
    def __init__(self, root_paths, temperature, volume, steps=None, **kwargs):
        """Parameters
        ---------
        root_paths: a list of strings
          Paths to ProtoMS output directories. Each string specifies an
          the output directory containing a different repeat of a
          calculation. Reported results are averaged over the given repeats.
        temperature: float
          Simulation temperature in degrees Kelvin
        volume: float
          The volume of the simulated GCMC region
        steps: integer, optional
          Number of steps to use in fitting the GCI neural network. The default
          value of None attempts to guess the required number of steps.
        **kwargs:
            Additional keyword arguments are passed to Estimator classes
            during initialisation.
        """
        self.subdir = ''
        feb.FreeEnergyCalculation.__init__(
            self,
            root_paths=[root_paths],
            temperature=temperature,
            estimators=[GCI],
            volume=volume,
            steps=steps,
            **kwargs)

    def _path_constructor(self, root_path):
        """Given a root_path (string) construct a full path
        suitable for globbing to find output directories."""
        return os.path.join(root_path, "b_*", self.subdir)

    def _get_lambda(self, path):
        """Given a path (string) extract the contained lambda value"""
        return float(path.split('/')[-2].split('_')[1])

    def calculate(self, subset=(0., 1., 1)):
        """For each estimator return the evaluated potential of mean force.

        Parameters
        ----------
        subset: tuple of two floats and an int
            specify the subset of data to use in the calculation
        """
        results = {}
        for est, legs in self.estimators.items():
            leg_result = GCMCResult()
            for i, leg in enumerate(legs):
                leg_result += GCMCResult([
                    rep.subset(*subset).calculate(self.temperature)
                    for rep in leg
                ])
            results[est] = leg_result
        return results[GCI]

    def _body(self, args):
        """The business logic of the calculation.

        Parameters
        ----------
        args: argparse.Namespace object
            Namespace from argumentparser
        """
        results = self.calculate(subset=(args.lower_bound, args.upper_bound))

        fig, ax = plt.subplots()
        plot_titration(results, ax)
        self.figures['titration'] = fig

        fig, table = plot_insertion_pmf(results)
        self.tables.append(table)
        self.figures['insertion_pmf'] = fig

        eqb_B = results.equilibrium_B
        eqb_index = np.argmin(
            [abs(B - eqb_B) for B in results.occupancies.coordinate])
        closest_B = results.occupancies.coordinate[eqb_index]

        self.footer += 'The equilibrium B value is %.3f\n' % eqb_B
        self.footer += 'Most similar simulated B value is %.3f\n' % closest_B
        self.footer += 'Occupancy at %.3f is %s\n' % \
                       (closest_B, results.occupancies.values[eqb_index])

        min_index = np.argmin(results.insertion_pmf.values)
        self.footer += "\nOccupancy at insertion PMF minimum is %.3f\n" % \
                       results.insertion_pmf.coordinate[min_index]
        return results


def plot_titration(results, ax, dot_fmt='b'):
    """Convenience function to plot the titration data from repeats."""
    for rep in results.data[0]:
        ax.plot(rep.coordinate, rep.values, 'o')
    results.model.plot(ax, xlabel='B Value', ylabel='Occupancy', color='black')


def plot_insertion_pmf(results, title=''):
    """Convenience function to plot the insertion pmf data from repeats"""
    table = Table(
        title,
        fmts=['%d', '%.3f', '%.3f', '%.3f'],
        headers=['Number of Waters',
                 'Insertion Free Energy',
                 'Network Binding Free Energy',
                 'Water Binding Free Energy'])

    steps = results.data[0][0].pmf.coordinate
    pmf = feb.PMF(steps, *[rep.pmf for rep in results.data[0]])
    prev_fe = 0.0
    for i, fe in enumerate(pmf):
        bind_fe = fe - i*tip4p_excess
        table.add_row([i, fe, bind_fe, bind_fe - prev_fe])
        prev_fe = bind_fe.value

    fig, ax = plt.subplots()
    pmf.plot(ax, xlabel="Occupancy")
    return fig, table


def get_arg_parser():
    """Add custom argument parser for this script"""
    parser = feb.FEArgumentParser(
        description="Calculate water binding free energies using Grand "
                    "Canonical Integration.",
        parents=[feb.get_base_arg_parser()],
        conflict_handler='resolve')
    parser.add_argument(
        '-d', '--directories', nargs='+', required=True,
        help="Location of folders containing ProtoMS output subdirectories. "
             "Multiple directories can be supplied to this flag and indicate "
        "repeats of the same calculation.")
    parser.add_argument(
        '-v', '--volume', required=True, type=float,
        help="Volume of the calculations GCMC region.")
    parser.add_argument(
        '-n', '--nsteps', type=int,
        help='Override automatic guessing of maximum number of waters for '
             'titration curve fitting.')
    return parser


def run_script(cmdline):
    """Execute the script, allows for straight-forward testing."""
    args = get_arg_parser().parse_args(cmdline)
    tc = TitrationCalculation(
        args.directories,
        args.temperature,
        args.volume,
        args.nsteps,
        results_name=args.name)
    tc.run(args)


if __name__ == '__main__':
    run_script(sys.argv[1:])
