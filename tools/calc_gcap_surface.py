from glob import glob
import numpy as np
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import os
from scipy import optimize
import sys
from protomslib import free_energy as fe

if "DISPLAY" not in os.environ or os.environ["DISPLAY"] == "":
    matplotlib.use("Agg")
import matplotlib.pyplot as plt


class SurfaceCalculation(fe.FreeEnergyCalculation):
    """Calculate the free energy surface arising from a two-dimensional
    GCAP simulation.
    """

    def __init__(
        self,
        root_path,
        temperature,
        volume,
        results_name="results",
        estimators=[fe.TI, fe.BAR, fe.MBAR],
    ):
        """Parameters
        ---------
        root_path: string
          The path of the root directory containing GCAP simulation data
        temperature: float
          The simulated temperature in Kelvin
        volume: float
          Volume of the GCMC region in Angstrom^3
        results_name: string
          Filename of the ProtoMS results file to use
        estimators:  list of estimator classes, optional
          The estimators to use to calculate pmfs across lambda
        """
        fe.FreeEnergyCalculation.__init__(
            self, [root_path], temperature, extract_data=False, volume=volume
        )
        self.paths = self.paths[0][0]

        # determine simulation B values and store them
        self.B_values = list(
            map(self._get_bvalue, glob("%s/*" % self.paths[0]))
        )
        self.B_values.sort()

        self.gci_calculations = [
            fe.TitrationCalculation(
                ["%s/lam-%.3f" % (root_path[0], lam)],
                temperature,
                volume,
                results_name=results_name,
                nfits=5,
            )
            for lam in self.lambdas[0][0]
        ]

        self.fep_calculations = [
            fe.FreeEnergyCalculation(
                [root_path],
                temperature,
                subdir="b_%.3f" % B,
                estimators=estimators,
                results_name=results_name,
            )
            for B in self.B_values
        ]

    def _get_bvalue(self, path):
        """Given a path (string) extract the contained B value"""
        return float(path.split("/")[-1].split("_")[1])

    def _body(self, args):
        """Calculation business logic.

        Parameters
        ----------
        args: argparse.Namespace object
            Namespace from argumentparser
        """
        fep_results = [calc.calculate() for calc in self.fep_calculations]
        gci_results = [calc.calculate() for calc in self.gci_calculations]

        # calculate pmfs across lambda
        gci_pmfs = []
        for i, result in enumerate(gci_results):
            model = result.data[0][0].model
            pmf = []
            N_min = min(result.occupancies.as_floats())
            for n in result.occupancies.as_floats():
                trans_nrg = fe.insertion_pmf(
                    np.array([N_min, n]), model, args.volume
                )[1]
                pmf.append(trans_nrg - n * fe.tip4p_excess)
            gci_pmfs.append(pmf)
        gci_pmfs = np.array(gci_pmfs)

        # for each alchemical estimator calculate pmfs across B
        # and construct the free eneergy surface
        results = {}
        lams, Bs = self.lambdas[0][0], self.B_values
        for est in self.estimators:
            fep_pmfs = np.array(
                [result[est].pmf.as_floats() for result in fep_results]
            ).T

            surface = np.zeros(len(self.lambdas[0][0]) * len(self.B_values))
            surf, fopt = optimize.fmin_bfgs(
                _obj_func,
                surface,
                args=(gci_pmfs, fep_pmfs),
                disp=False,
                full_output=True,
            )[:2]

            surf = surf.reshape(fep_pmfs.shape)
            fig = plot_surface(
                lams,
                Bs,
                surf,
                "dA (kcal/mol)",
                "Free Energy Surface - %s" % est.__name__,
            )
            self.figures["pmf2d_%s" % est.__name__] = fig
            results[est] = surf
        fig = plot_surface(
            lams,
            Bs,
            np.array(
                [result.occupancies.as_floats() for result in gci_results]
            ),
            "Number of Waters",
            "Water Occupancy Surface",
        )
        self.figures["occupancy2d"] = fig
        return results


def _obj_func(x, y, z):
    """Objective function minimised to produce free energy surface
    x is the trial surface, y is the free energy slices over B and z are
    the free energy slices over lambda.

    The objective function compares relative heights of neighbouring points
    on the trial surface with those of the input data free energy curves.
    """
    x = x.reshape(y.shape)
    tot = 0
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            if i == 0 and j == 0:
                pass
            elif i == 0:
                tot += ((x[i, j] - x[i, j - 1]) - (y[i, j] - y[i, j - 1])) ** 2
            elif j == 0:
                tot += ((x[i, j] - x[i - 1, j]) - (z[i, j] - z[i - 1, j])) ** 2
            else:
                tot += (
                    (x[i, j] - x[i - 1, j]) - (z[i, j] - z[i - 1, j])
                ) ** 2 + (
                    (x[i, j] - x[i, j - 1]) - (y[i, j] - y[i, j - 1])
                ) ** 2
    return tot


def plot_surface(lambdas, Bs, Z, zlabel="", title=""):
    """Utility function to plot a 2D free energy surface.

    Parameters
    ----------
    lambdas: list-like of numbers
      the simulated lambda values
    Bs: list-like of numbers
      the simulated B values
    Z: 2D list-like of numbers
      the heights of the free energy surface
    """
    fig = plt.figure()
    ax = Axes3D(fig)
    Bs, lambdas = np.meshgrid(Bs, lambdas)
    ax.plot_surface(Bs, lambdas, Z, cstride=1, rstride=1)

    ax.set_xlabel("B")
    ax.set_ylabel("lambda")
    ax.set_zlabel(zlabel)
    ax.set_title(title)

    return fig


def get_arg_parser():
    """Returns the custom argument parser for this script"""
    parser = fe.FEArgumentParser(
        description="Calculate free energy differences using a range of"
        " estimators",
        parents=[fe.get_alchemical_arg_parser()],
    )
    parser.add_argument(
        "-v",
        "--volume",
        type=float,
        default=None,
        help="Volume of GCMC region",
    )
    parser.add_argument(
        "--estimators",
        nargs="+",
        default=["ti", "mbar", "bar"],
        choices=["ti", "mbar", "bar"],
        help="Choose free energy estimator to use. By default TI, BAR and MBAR"
        " are used. Note that the GCAP estimator assumes a different file"
        " structure and ignores the --subdir flag.",
    )
    return parser


def run_script(cmdline):
    """Execute the script, allows for straight-forward testing."""
    class_map = {"ti": fe.TI, "mbar": fe.MBAR, "bar": fe.BAR}
    args = get_arg_parser().parse_args(cmdline)
    calc = SurfaceCalculation(
        args.directories[0],
        args.temperature,
        estimators=list(map(class_map.get, args.estimators)),
        volume=args.volume,
        results_name=args.name,
    )
    calc.run(args)


if __name__ == "__main__":
    run_script(sys.argv[1:])
