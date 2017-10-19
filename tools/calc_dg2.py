import free_energy_base


def get_arg_parser():
    """Add custom options for this script"""
    parser = free_energy_base.get_arg_parser()
    return parser


if __name__ == '__main__':
    args = get_arg_parser().parse_args()
    calc = free_energy_base.FreeEnergyCalculation(root_paths=args.directories)
    results = calc.calculate()

    for estimator in results:
        print estimator
        for pmf, path in zip(results[estimator], args.directories):
            print "%s: %.4f" % (path, pmf.dG)
        print
