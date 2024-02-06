import matplotlib
import numpy as np
import os
import pymbar
from protomslib import free_energy as fe
from protomslib import simulationobjects as sim

if "DISPLAY" not in os.environ or os.environ["DISPLAY"] == "":
    matplotlib.use("Agg")
import matplotlib.pyplot as plt


class GCSEnergies(fe.Estimator):
    def __init__(self, lambdas, subdir_glob, **kwargs):
        fe.Estimator.__init__(
            self,
            lambdas,
            subdir_glob=subdir_glob,
            results_name="results_inst",
            **kwargs
        )

    def add_data(self, series):
        tmp = []
        for term, dat in series.interaction_energies.items():
            if "GCS" in term:
                tmp.append(dat[0].curr)
        self.data.append(np.sum(tmp, axis=0))
        return len(self.data[-1])


class GCSLongEnergies(GCSEnergies):
    # results_name = 'results_long'
    def __init__(self, lambdas, subdir_glob, **kwargs):
        fe.Estimator.__init__(
            self,
            lambdas,
            subdir_glob=subdir_glob,
            results_name="results_long",
            **kwargs
        )


class ReweightedTitration(fe.TitrationCalculation):
    def __init__(
        self,
        root_paths,
        temperature,
        volume,
        steps,
        filename,
        longname,
        subdir="",
    ):
        fe.FreeEnergyCalculation.__init__(
            self,
            root_paths=root_paths,
            temperature=temperature,
            estimators=[fe.GCI, GCSEnergies, GCSLongEnergies],
            subdir=subdir,
            extract_data=False,
        )
        self.volume = volume
        self.steps = steps
        # adjust estimator instance attributes so that the
        # correct filenames are read, this is a bit ugly
        for est in self.estimators[fe.GCI][0]:
            est.results_name = filename
        for est in self.estimators[GCSEnergies][0]:
            est.results_name = filename
        for est in self.estimators[GCSLongEnergies][0]:
            est.results_name = longname
        self._extract_series(self.paths)

    def calculate(self, subset=(0.0, 1.0, 1)):
        GCI_estimators = {"raw": [], "umbrella": [], "mbar": []}
        beta = 1 / (sim.boltz * self.temperature)

        for i, rep in enumerate(self.estimators[fe.GCI][0]):
            # get relevant quantities into nice numpy arrays
            # occupancy data
            ns = np.array(rep.subset(*subset).data)
            # short cutoff interaction energies
            es = np.array(
                self.estimators[GCSEnergies][0][i].subset(*subset).data
            )
            # long cutoff interaction energies
            les = np.array(
                self.estimators[GCSLongEnergies][0][i].subset(*subset).data
            )
            # simulated B values
            Bs = np.array(map(self._get_lambda, self.paths[0][i]))
            # simulated chemical potentials
            mus = [
                B_to_chemical_potential(B, self.volume, self.temperature)
                for B in Bs
            ]

            # put quantities together to make u_kn, see MBAR documentation
            N = ns.shape[1]
            Nsims = len(mus)
            u_kn = np.zeros(shape=(2 * Nsims, N * Nsims))
            for i in range(Nsims):
                for j in range(Nsims):
                    u_kn[i, j * N : (j + 1) * N] = beta * (
                        es[j] + ns[j] * mus[i]
                    )
                    u_kn[i + Nsims, j * N : (j + 1) * N] = beta * (
                        les[j] + ns[j] * mus[i]
                    )

            # Create array of sample numbers. Entries with long-range
            # energies are left as 0
            N_k = np.zeros(Nsims * 2, dtype=int)
            N_k[:Nsims] = N

            # MBAR reweighting
            mbar = pymbar.MBAR(u_kn, N_k)
            mbar_ns = mbar.computeExpectations(ns.reshape(Nsims * N))[0]

            # also do reweighting with conventional umbrella result
            umbrella_ns = np.zeros_like(mus)
            for i in range(Nsims):
                exps = np.exp(beta * (es[i] - les[i]))
                umbrella_ns[i] = (exps * ns[i]).mean() / exps.mean()

            # rep is a suitable GCI estimator for raw occupancy data
            GCI_estimators["raw"].append(rep.subset(*subset))

            # create GCI estimators to calculate binding free energies
            # from reweighted occupancy data
            # for some reason reweighted data comes out of MBAR in reverse
            mbar = fe.GCI(Bs)
            mbar.data = mbar_ns[Nsims:][::-1].reshape((Nsims, 1))
            GCI_estimators["mbar"].append(mbar)

            umbrella = fe.GCI(Bs)
            umbrella.data = umbrella_ns.reshape((Nsims, 1))
            GCI_estimators["umbrella"].append(umbrella)

        # gather calculatinos from each repeat and package into
        # GCIResult objects for different reweighting methods.
        results = {}
        for key in "raw", "mbar", "umbrella":
            results[key] = fe.GCMCResult(
                [
                    rep.calculate(self.temperature, self.volume, self.steps)
                    for rep in GCI_estimators[key]
                ]
            )
        return results

    def _body(self, args):
        results = self.calculate(subset=(args.lower_bound, args.upper_bound))
        keys = ("raw", "mbar", "umbrella")

        # plot titration curves and fitted models
        self.figures["titration"], ax = plt.subplots()
        for key in keys:
            line = results[key].occupancies.plot(ax, fmt="o", show_error=False)
            results[key].model.plot(
                ax,
                color=line.get_color(),
                label=key,
                xlabel="B Value",
                ylabel="Occupancy",
            )
        ax.legend(loc="best")

        # plot insertion_pmfs
        self.figures["insertion_pmf"], ax = plt.subplots()
        for key in keys:
            results[key].insertion_pmf.plot(ax, label=key, xlabel="Occupancy")
        ax.legend(loc="best")

        # print out binding free energies in table format
        for key in keys:
            table = fe.Table(
                key.upper(),
                fmts=["%d", "%.3f"],
                headers=["Number of Waters", "Binding Free Energy"],
            )
            for i, dA in enumerate(results[key].insertion_pmf):
                table.add_row([i, dA])
            self.tables.append(table)
        return results


def B_to_chemical_potential(B, volume, temperature):
    beta = 1 / (temperature * sim.boltz)
    thermal_wavelength = 1.00778365325  # of water, in angstroms
    return (B - np.log(volume / thermal_wavelength**3)) / beta


def get_arg_parser():
    """Add custom options for this script"""
    parser = fe.FEArgumentParser(
        description="Calculate water binding free energies using Grand "
        "Canonical Integration and incorporating rigorous "
        "corrections for long-range electrostatics and "
        "dispersion. Reweighting is carried out using both "
        "an umbrella sampling result and MBAR.",
        parents=[fe.get_gci_arg_parser()],
        conflict_handler="resolve",
    )
    parser.add_argument(
        "-f",
        "--filename",
        default="results_inst",
        help="Name of file containing simulation energies with short cutoff."
        " Default=results_inst",
    )
    parser.add_argument(
        "-ln",
        "--longname",
        default="results_long",
        help="Name of file containing simulation energies with long cutoff."
        " Default=results_long",
    )
    return parser


if __name__ == "__main__":
    args = get_arg_parser().parse_args()
    tc = ReweightedTitration(
        [args.directories],
        args.temperature,
        args.volume,
        args.nsteps,
        args.filename,
        args.longname,
        subdir=args.subdir,
    )
    tc.run(args)
    plt.show()
