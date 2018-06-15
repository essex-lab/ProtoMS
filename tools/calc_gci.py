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
        parents=[fe.get_gci_arg_parser()],
        conflict_handler='resolve')
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
