import matplotlib
import numpy as np
import os
import free_energy_base as feb

if "DISPLAY" not in os.environ or os.environ["DISPLAY"] == "":
    matplotlib.use('Agg')
import matplotlib.pyplot as plt


class FreeEnergyCalculation(feb.FreeEnergyCalculation):
    def test_equilibration(self, discard_limit):
        """Perform a series of calculations discarding increasing portions
        from the start of the loaded data series. This assesses whether
        unequilibrated data from the start of a simulation is biasing
        the result.

        Parameters
        ----------
        discard_limit : float
          Value between 0 and 1 that indicates the maximum amount of
          data to discard. E.g. a value of 0.2 would perform a series
          of calculations that range between using all of available
          data and the last 80% of the data.

        Returns
        -------
        dict
          keys: data proportions
          values: outputs from self.calculate(subset=(X, 1.0))
                  where X is the corresponding key
        """
        proportions = np.linspace(0.0, discard_limit, 10)
        return {prop: self.calculate(subset=(prop, 1.0))
                for prop in proportions}

    def test_convergence(self, discard_limit, lower_limit=0.):
        """Perform a series of calculations discarding increasing portions
        from the end of the loaded data series. This assesses whether
        sufficient data has been sampled to produce converged free energy
        estimates.

        Parameters
        ----------
        discard_limit : float
          Value between 0 and 1 that indicates the minimum proportion of data
          to use. E.g. a value of 0.8 would perform a series of calculations
          that range between using the all of the available data and the first
          80% of the data.
        lower_limit : float, optional
          Specifies an initial proportion of the data set to discard from all
          calculations. Allows unequilibrated data to be removed.

        Returns
        -------
        dict
          keys: data proportions
          values: outputs from self.calculate(subset=(lower_limit, X))
                  where X is the corresponding key
        """
        proportions = np.linspace(discard_limit, 1.0, 10)
        return {prop: self.calculate(subset=(lower_limit, prop))
                for prop in proportions}

    def _body(self, args):
        """Execute the required calculation according to the values present in
        the argparse.Namespace args.
        """
        if (args.test_equilibration or args.test_convergence) is not None:
            if (args.test_equilibration) is not None:
                results = self.test_equilibration(args.test_equilibration)
            else:
                results = self.test_convergence(args.test_convergence,
                                                lower_limit=args.lower_bound)
        else:
            results = self.calculate(subset=(args.lower_bound,
                                             args.upper_bound))

        if (args.test_equilibration or args.test_convergence) is not None:
            fig = plot_fractional_dataset_results(results, calc.estimators)
            figname = 'equil' if args.test_equilibration else 'convergence'
            self.figures[figname] = fig
        else:
            if args.pmf:
                self.figures['pmf'] = plot_pmfs(results)
            # print_results(args.directories, results)
            self.tables.extend(results_tables(args.directories, results))
        return results


def plot_results(x, results, ax, **kwargs):
    y = np.array([result.dG.value for result in results])
    err = np.array([result.dG.error for result in results])

    line = ax.plot(x, y, **kwargs)[0]
    ax.plot(x, y+err, '--', linewidth=1, color=line.get_color())
    ax.plot(x, y-err, '--', linewidth=1, color=line.get_color())


def plot_fractional_dataset_results(results, estimators):
    """Graph results of calculations that use variable portions of available
    data i.e. equilibration and convergence test.s."""
    fig, ax = plt.subplots()
    for estimator in estimators:
        x = sorted(results)
        dat = [results[prop][estimator] for prop in x]
        plot_results(x, dat, ax, label=estimator.__name__)
    ax.legend(loc='best')
    ax.set_xlabel('Proportion')
    ax.set_ylabel('Free energy (kcal/mol)')
    return fig


def plot_pmfs(results):
    """Graph average potentials of mean force for all estimators."""
    fig, ax = plt.subplots()
    for estimator in sorted(results):
        results[estimator].pmf.plot(ax, label=estimator.__name__)
    ax.legend(loc='best')
    ax.set_xlabel('Lambda value')
    ax.set_ylabel('Free energy (kcal/mol)')
    return fig


# def print_results(directories, results):
#     """Print calculated free energies. If multiple repeats are present
#     the mean and standard error are also printed.
#     """
#     for estimator in sorted(results, key=lambda x: x.__name__):
#         print "%s -" % estimator.__name__
#         for i, res in enumerate(results[estimator].data):
#             dGs = []
#             for j, pmf in enumerate(res):
#                 dGs.append(pmf.dG)
#                 print "%-20s: %.4f" % (directories[i][j], pmf.dG)
#             print "%-20s: %s\n" % ('Mean', feb.FreeEnergy.fromData(dGs))
#         print "Total Mean: %s\n" % results[estimator].dG

def results_tables(directories, results):
    """Print calculated free energies. If multiple repeats are present
    the mean and standard error are also printed.
    """
    tables = []
    for estimator in sorted(results, key=lambda x: x.__name__):
        table = feb.Table(estimator.__name__, ['%s', "%.2f"])
        # print "%s -" % estimator.__name__
        for i, res in enumerate(results[estimator].data):
            dGs = []
            for j, pmf in enumerate(res):
                dGs.append(pmf.dG)
                table.add_row([directories[i][j], pmf.dG])
            table.add_row(['Mean', feb.FreeEnergy.fromData(dGs)])
            table.add_blank_row()
        table.add_row(['Total Mean', results[estimator].dG])
        tables.append(table)
    return tables


def get_arg_parser():
    """Add custom options for this script"""
    parser = feb.FEArgumentParser(
        description="Calculate free energy differences using a range of"
                    " estimators",
        parents=[feb.get_arg_parser()])
    parser.add_argument(
        '--pmf', action='store_true', default=False,
        help="Make graph of potential of mean force",
        clashes=('test_convergence', 'test_equilibration'))
    parser.add_argument(
        '--test-equilibration', default=None, type=float,
        help="Perform free energy calculations 10 times using varying "
             "proportions of the total data set provided. Data used will "
             "range from 100%% of the dataset down to the proportion "
             "provided to this argument",
        clashes=('test_convergence', 'lower_bound'))
    parser.add_argument(
        '--test-convergence', default=None, type=float,
        help="Perform free energy calculations 10 times using varying "
             "proportions of the total data set provided. Data used will "
             "range from 100%% of the dataset up to the proportion "
             "provided to this argument")
    return parser


if __name__ == '__main__':
    args = get_arg_parser().parse_args()
    calc = FreeEnergyCalculation(root_paths=args.directories,
                                 temperature=args.temperature,
                                 subdir=args.subdir)
    calc.run(args)
