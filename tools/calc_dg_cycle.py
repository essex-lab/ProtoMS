from collections import defaultdict
from glob import glob
import numpy as np
import free_energy_base as feb


class FreeEnergy(object):
    def __init__(self, value, error):
        self.value = value
        self.error = error

    def __add__(self, other):
        return FreeEnergy(self.value + other.value,
                          (self.error**2 + other.error**2)**0.5)

    def __sub__(self, other):
        return FreeEnergy(self.value - other.value,
                          (self.error**2 + other.error**2)**0.5)

    def __str__(self):
        return "%.4f +/- %.4f" % (self.value, self.error)

    def __repr__(self):
        return "<FreeEnergy: value=%.4f error=%.4f>" % (self.value, self.error)


def get_arg_parser():
    """Add custom options for this script"""
    parser = feb.FEArgumentParser(
        description="Calculate free energy differences using a range of"
                    " estimators",
        parents=[feb.get_arg_parser()])
    parser.add_argument(
        "-s", "--signs", nargs='+', type=str, choices=('+', '-'),
        help="", required=True)

    # use mutually exclusive group for these as we need one and only one
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--dualtopology", action='store_true', default=False,
        help="")
    group.add_argument(
        "--singletopology", nargs=1, type=str, choices=('comp', 'decomp'),
        help="")
    return parser


def print_solv_bind(dG_solvs, dG_binds, directories, signs, estimator):
    fmt = "%20s  %20s  %20s"
    print fmt % ('', 'ddG Solvation', 'ddG Binding')
    closure_solv = FreeEnergy(0., 0.)
    closure_bind = FreeEnergy(0., 0.)
    for root, sign in zip(directories, signs):
        dG_solv = dG_solvs[root][estimator]
        dG_bind = dG_binds[root][estimator]
        print fmt % (root, dG_solv, dG_bind)
        if sign == '+':
            closure_solv += dG_solv
            closure_bind += dG_bind
        else:
            closure_solv -= dG_solv
            closure_bind -= dG_bind
    print fmt % ('Cycle Closure', closure_solv, closure_bind)
    print


if __name__ == '__main__':
    args = get_arg_parser().parse_args()

    if len(args.directories) != len(args.signs):
        raise Exception(
            "Please provide one valid sign for each provided directory.")

    output_dir_format = "%s/out?%s_%s"

    if args.singletopology == 'decomp':
        midfixes = ['_ele', '_vdw']
    elif args.singletopology == 'comp':
        midfixes = ['_comp']
    else:
        midfixes = ['']

    legs = ('gas', 'free', 'bnd')
    estimators = (feb.TI, feb.BAR, feb.MBAR)
    results = {}

    for root in args.directories:
        results[root] = {}
        for leg in legs:
            leg_dGs = {est: FreeEnergy(0., 0.) for est in estimators}
            for mid in midfixes:
                output_dir = output_dir_format % (root, mid, leg)
                calc = feb.FreeEnergyCalculation(
                    root_paths=glob(output_dir),
                    temperature=args.temperature,
                    estimators=estimators)
                data = calc.calculate(
                    subset=(args.lower_bound, args.upper_bound, 1))
                for est in estimators:
                    dat = np.array([pmf.dG for pmf in data[est]])
                    if len(dat) == 0:
                        fe = FreeEnergy(float('nan'), float('nan'))
                    else:
                        fe = FreeEnergy(dat.mean(), dat.std()/len(dat)**0.5)
                    leg_dGs[est] += fe
            results[root][leg] = leg_dGs

    dG_solvation = {}
    dG_binding = {}
    for root in args.directories:
        dG_solvation[root] = {}
        dG_binding[root] = {}
        for est in estimators:
            gas = results[root]['gas'][est]
            free = results[root]['free'][est]
            bnd = results[root]['bnd'][est]
            if args.dualtopology:
                dG_solvation[root][est] = free
                dG_binding[root][est] = bnd - free
            else:
                dG_solvation[root][est] = free - gas
                dG_binding[root][est] = bnd - free

    for est in estimators:
        print "%s -" % est.__name__
        fmt = "%20s  %20s  %20s  %20s"
        print fmt % ('', 'dG gas', 'dG free', 'dG bound')
        # closures = defaultdict(lambda: FreeEnergy(0., 0.))
        closures = {leg: FreeEnergy(0., 0.) for leg in legs}
        for root, sign in zip(args.directories, args.signs):
            sub_vals = (root,)
            for leg in legs:
                dG = results[root][leg][est]
                sub_vals += (dG,)
                if sign == '+':
                    closures[leg] += dG
                else:
                    closures[leg] -= dG
            print fmt % sub_vals
        print fmt % (('Cycle Closure',) + tuple(closures[leg] for leg in legs))
        print

        print_solv_bind(dG_solvation, dG_binding,
                        args.directories, args.signs, est)
