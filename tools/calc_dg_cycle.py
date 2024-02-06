from glob import glob
import sys
from protomslib import free_energy as fe


def state_data_table(results, directories, signs, states, estimator):
    """Returns a Table object containing the free energy differences for the
    individual states i.e. bound, free and gas.

    Parameters
    ----------
    results: dict such that results[root][state][estimator] = Result object
      The free energy results to tabulate
    directories: list of strings
      List of root directories containing calculation data
    signs: list of '+' or '-'
      The signs to use when calculating the free energy cycle
    states: list of strings
      The thermodynamic states simulated
    estimator: Estimator Class
      The estimator to tabulate results of
    """
    table = fe.Table(
        estimator.__name__,
        fmts=["%s:", "%.3f", "%.3f", "%.3f"],
        headers=["", "dG gas", "dG free", "dG bound"],
    )

    closures = {state: fe.Quantity(0.0, 0.0) for state in states}
    for root, sign in zip(directories, signs):
        root_dGs = (root,)
        for state in states:
            dG = results[root][state][estimator]
            root_dGs += (dG,)
            if sign == "+":
                closures[state] += dG
            else:
                closures[state] -= dG
        table.add_row(root_dGs)
    table.add_row(["Cycle Closure"] + [closures[state] for state in states])
    return table


def solv_bind_table(dG_solvs, dG_binds, directories, signs, estimator):
    """Returns a Table object containing free energies of solvation and
    binding.

    Parameters
    ----------
    dG_solvs: dictionary[root][estimator] = Result object
      Calculated solvation free energy results
    dG_binds: dictionary[root][estimator] = Result object
      Calculated binding free energy results
    directories: list of strings
      List of root directories containing calculation data
    signs: list of '+' or '-'
      The signs to use when calculating the free energy cycle
    estimator: Estimator Class
      The estimator to tabulate results of
    """
    table = fe.Table(
        "",
        fmts=["%s:", "%.3f", "%.3f"],
        headers=["", "ddG Solvation", "ddG Binding"],
    )

    closure_solv = fe.Quantity(0.0, 0.0)
    closure_bind = fe.Quantity(0.0, 0.0)
    for root, sign in zip(directories, signs):
        dG_solv = dG_solvs[root][estimator]
        dG_bind = dG_binds[root][estimator]
        table.add_row([root, dG_solv, dG_bind])
        if sign == "+":
            closure_solv += dG_solv
            closure_bind += dG_bind
        else:
            closure_solv -= dG_solv
            closure_bind -= dG_bind
    table.add_row(["Cycle Closure", closure_solv, closure_bind])
    return table


class CycleCalculation(fe.FreeEnergyCalculation):
    """Perform multiple free energy calculations and combine the
    results to calculate a cycle.
    """

    def __init__(
        self,
        estimators=[fe.TI, fe.BAR, fe.MBAR],
        volume=None,
        results_name="results",
    ):
        """Parameters
        ---------
        estimators: list of estimator classes, optional
          The estimators to use
        volume: float, optional
          Volume of the GCMC region
        results_name: string, optional
          Filename of the ProtoMS results file to use
        """
        self.estimators = estimators
        self.figures = {}
        self.tables = []
        self.volume = volume
        self.results_name = results_name
        self.header = ""
        self.footer = ""

    def _body(self, args):
        """Calculation business logic.

        Parameters
        ----------
        args: argparse.Namespace object
            Namespace from argumentparser
        """
        if len(args.directories) != len(args.signs):
            raise Exception(
                "Please give exactly one sign for each provided directory."
            )

        # determine what output folder names we're looking for
        # elements of midfixes will replace the central %s of
        # the output_dir_format string
        output_dir_format = "%s/out?%s_%s"
        if args.singletopology == "sep":
            # calculation is split into electrostatics and van der Waals
            # components, need to do both and add them up
            midfixes = ["_ele", "_vdw"]
        elif args.singletopology == "comb":
            midfixes = ["_comb"]
        else:
            midfixes = [""]

        # perform calculations
        # end up with results dictionary structured as -
        # results[root][state][estimator] == fe.Quantity object
        # thus we have one free energy difference for each root
        # directory, for each state (gas, free or bound), calculated
        # with each of the estimators
        states = ("gas", "free", "bnd")
        results = {}
        for root in args.directories:
            results[root] = {}
            for state in states:
                state_dGs = {
                    est: fe.Quantity(0.0, 0.0) for est in self.estimators
                }
                for mid in midfixes:
                    output_dir = output_dir_format % (root, mid, state)
                    calc = fe.FreeEnergyCalculation(
                        root_paths=[glob(output_dir)],
                        temperature=args.temperature,
                        subdir=args.subdir,
                        estimators=self.estimators,
                        volume=self.volume,
                        results_name=self.results_name,
                    )
                    data = calc.calculate(
                        subset=(args.lower_bound, args.upper_bound, 1)
                    )
                    for est in self.estimators:
                        state_dGs[est] += data[est].dG

                results[root][state] = state_dGs

        # combine state changes to get meaningful free energy differences
        dG_solvation, dG_binding = {}, {}
        for root in args.directories:
            dG_solvation[root] = {}
            dG_binding[root] = {}
            for est in self.estimators:
                gas = results[root]["gas"][est]
                free = results[root]["free"][est]
                bnd = results[root]["bnd"][est]

                # solvation depends on calculation type
                if args.dualtopology:
                    dG_solvation[root][est] = free
                else:
                    dG_solvation[root][est] = free - gas
                dG_binding[root][est] = bnd - free

        for est in self.estimators:
            self.tables.append(
                state_data_table(
                    results, args.directories, args.signs, states, est
                )
            )
            # print
            self.tables.append(
                solv_bind_table(
                    dG_solvation, dG_binding, args.directories, args.signs, est
                )
            )

        return results


def get_arg_parser():
    """Returns the custom argument parser for this script"""
    parser = fe.FEArgumentParser(
        description="High level script that attempts to use data from multiple"
        " calculations to provide free energies of solvation and "
        "binding. Also calculates cycle closures for all data. "
        "Assumes standard ProtoMS naming conventions for data "
        "output directories. Data should be organised such that "
        "each transformation between two ligands should have a "
        "single master directory containing output directories "
        "for each simulation state (e.g. master/out1_free). Master"
        "directories should be passed to the -d flag. "
        "Reported free energies are averages "
        "over all repeats found. Reported errors are single "
        "standard errors calculated from repeats.",
        parents=[fe.get_alchemical_arg_parser()],
        conflict_handler="resolve",
    )
    parser.add_argument(
        "-d",
        "--directories",
        nargs="+",
        required=True,
        help="Location of folders containing ProtoMS output directories.",
    )
    parser.add_argument(
        "-s",
        "--signs",
        nargs="+",
        type=str,
        choices=("+", "-"),
        help="List of '+' or '-' characters, one for each directory provided "
        "to the -d flag. Indicates the sign that should be used for each "
        "free energy difference when calculating cycle closures.",
        required=True,
    )
    parser.add_argument(
        "--estimators",
        nargs="+",
        default=["ti", "mbar", "bar"],
        choices=["ti", "mbar", "bar", "gcap"],
        help="Choose estimators",
    )
    parser.add_argument(
        "-v",
        "--volume",
        type=float,
        default=None,
        help="Volume of GCMC region",
    )

    # use mutually exclusive group for these as we need one and only one
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--dualtopology",
        action="store_true",
        default=False,
        help="Indicates data is for a dual topology calculation.",
    )
    group.add_argument(
        "--singletopology",
        type=str,
        choices=("comb", "sep"),
        help="Indicates data is for a single topology calculation. "
        "Option comb indicates a single step calculation. "
        "Option sep indicates separate steps for van der Waals "
        " and electrostatics components.",
    )
    return parser


def run_script(cmdline):
    """Execute the script, allows for straight-forward testing."""
    class_map = {
        "ti": fe.TI,
        "mbar": fe.MBAR,
        "bar": fe.BAR,
        "gcap": fe.GCMCMBAR,
    }
    args = get_arg_parser().parse_args(cmdline)
    calc = CycleCalculation(
        estimators=list(map(class_map.get, args.estimators)),
        volume=args.volume,
        results_name=args.name,
    )
    calc.run(args)


if __name__ == "__main__":
    run_script(sys.argv[1:])
