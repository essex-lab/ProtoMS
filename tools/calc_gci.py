"""
Routines to calculate free energies with GCMC as outlined in
J. Am. Chem. Soc., 2015, 137 (47), pp 14930-14943

Can be executed from the command line as a stand-alone program
"""

import sys
from protomslib import free_energy as fe


def get_arg_parser():
    """Add custom argument parser for this script"""
    parser = fe.FEArgumentParser(
        description="Calculate water binding free energies using Grand "
                    "Canonical Integration.",
        parents=[fe.get_base_arg_parser()],
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
        help='Override automatic guessing of the number of steps to fit for '
             'titration curve fitting.')
    parser.add_argument(
        '--nmin', type=int,
        help='Override automatic guessing of the minimum number of waters for '
             'tittration curve fitting.')
    parser.add_argument(
        '--nmax', type=int,
        help='Override automatic guessing of maximum number of waters for '
             'titration curve fitting.')
    parser.add_argument(
        '--nfits', type=int, default=10,
        help='The number of independent fitting attempts for the neural '
             'network occupancy model. Increasing the number of fits may '
             'help improve results for noisy data.')
    return parser


def run_script(cmdline):
    """Execute the script, allows for straight-forward testing."""
    args = get_arg_parser().parse_args(cmdline)
    tc = fe.TitrationCalculation(
        args.directories,
        args.temperature,
        args.volume,
        args.nsteps,
        args.nmin,
        args.nmax,
        args.nfits,
        results_name=args.name)
    tc.run(args)


if __name__ == '__main__':
    run_script(sys.argv[1:])
