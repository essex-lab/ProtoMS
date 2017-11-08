import matplotlib.pyplot as plt
import numpy as np
import free_energy_base


def get_arg_parser():
    """Add custom options for this script"""
    parser = free_energy_base.FEArgumentParser(
        description="Calculate free energy differences using a range of"
                    " estimators",
        parents=[free_energy_base.get_arg_parser()])
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
    return fig, ax


if __name__ == '__main__':
    args = get_arg_parser().parse_args()
    calc = free_energy_base.FreeEnergyCalculation(root_paths=args.directories)

    if args.autoeqb:
        calc.autoeqb()

    if (args.test_equilibration or args.test_convergence) is not None:
        if (args.test_equilibration) is not None:
            results = calc.test_equilibration(args.test_equilibration)
        else:
            results = calc.test_convergence(args.test_convergence,
                                            lower_limit=args.lower_bound)
        fig, ax = plot_fractional_dataset_results(results, calc.estimators)
        ax.set_xlabel('proportion')
        ax.set_ylabel('free energy (kcal/mol)')
    else:
        results = calc.calculate(subset=(args.lower_bound, args.upper_bound))
        for estimator in sorted(results):
            print estimator.__name__
            dGs = [pmf.dG for pmf in results[estimator]]
            for dG, path in zip(dGs, args.directories):
                print "%s: %.4f" % (path, dG)
            if len(dGs) > 1:
                print "Mean: %.4f +/- %.4f" % (np.mean(dGs), 
                                               np.std(dGs)/len(dGs)**0.5)
            print

    plt.show()
