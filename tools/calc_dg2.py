import matplotlib.pyplot as plt
import numpy as np
import free_energy_base


def get_arg_parser():
    """Add custom options for this script"""
    parser = free_energy_base.FEArgumentParser(
        description="Calculate free energy differences using a range of"
                    " estimators",
        parents=[free_energy_base.get_arg_parser()])
    parser.add_argument(
        '--pmf', action='store_true', default=False,
        help="Make graph of potential of mean force",
        clashes=('test_convergence', 'test_equilibration'))
    return parser


def plot_fractional_dataset_results(results, estimators):
    fig, ax = plt.subplots()
    ys = {}
    for estimator in estimators:
        ys[estimator] = []
        for prop in sorted(results):
            ys[estimator].append(
                np.mean([pmf.dG for pmf in results[prop][estimator]]))

    xs = sorted(results)
    for estimator in ys:
        ax.plot(xs, ys[estimator], label=str(estimator))
    ax.legend()
    ax.set_xlabel('Proportion')
    ax.set_ylabel('Free energy (kcal/mol)')
    return fig, ax


def plot_pmfs(results):
    fig, ax = plt.subplots()
    for estimator in sorted(results):
        lambdas = results[estimator][0].lambdas
        pmf_2d = np.array([pmf.values for pmf in results[estimator]])
        ax.plot(lambdas, pmf_2d.mean(axis=0), label=estimator.__name__)
    ax.legend()
    ax.set_xlabel('Lambda value')
    ax.set_ylabel('Free energy (kcal/mol)')
    return fig, ax


def print_results(results):
    for estimator in sorted(results):
        print estimator.__name__
        dGs = [pmf.dG for pmf in results[estimator]]
        for dG, path in zip(dGs, args.directories):
            print "%s: %.4f" % (path, dG)
        if len(dGs) > 1:
            print "Mean: %.4f +/- %.4f" % (np.mean(dGs),
                                           np.std(dGs)/len(dGs)**0.5)
        print


if __name__ == '__main__':
    args = get_arg_parser().parse_args()
    calc = free_energy_base.FreeEnergyCalculation(root_paths=args.directories)

    results = calc.run(args)

    if (args.test_equilibration or args.test_convergence) is not None:
        plot_fractional_dataset_results(results, calc.estimators)
    else:
        if args.pmf:
            plot_pmfs(results)
        print_results(results)

    plt.show()
