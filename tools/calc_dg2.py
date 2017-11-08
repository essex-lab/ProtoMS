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
    return fig, ax


if __name__ == '__main__':
    args = get_arg_parser().parse_args()
    calc = free_energy_base.FreeEnergyCalculation(root_paths=args.directories)

    # if args.autoeqb:
    #     calc.autoeqb()

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
        if args.pmf:
            fig, ax = plt.subplots()
            ax.set_xlabel('Lambda value')
            ax.set_ylabel('Free energy (kcal/mol)')

        results = calc.calculate(subset=(args.lower_bound, args.upper_bound))
        for estimator in sorted(results):
            print estimator.__name__
            dGs = [pmf.dG for pmf in results[estimator]]
            for dG, path in zip(dGs, args.directories):
                print "%s: %.4f" % (path, dG)
            if args.pmf:
                lambdas = results[estimator][0].lambdas
                pmf_2d = np.array([pmf.values for pmf in results[estimator]])
                ax.plot(lambdas, pmf_2d.mean(axis=0), label=estimator.__name__)

                # print "\nPMF -"
                # print "Lambda   Free Energy (kcal/mol) "
                # for lam, val in zip(lambdas, pmf_2d.T):
                #     print "%.3f    %.4f" % (lam, val.mean())
            if len(dGs) > 1:
                print "Mean: %.4f +/- %.4f" % (np.mean(dGs),
                                               np.std(dGs)/len(dGs)**0.5)
            print

    if args.pmf:
        ax.legend()
            
    plt.show()
