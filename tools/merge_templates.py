# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
"""
Routine to merge a series of ProtoMS template files

This module defines a single public function:
merge_templates

Can be executed from the command line as a stand-alone program
"""

import logging
from protomslib import simulationobjects
from protomslib.templates import merge_templates

logger = logging.getLogger('protoms')


def get_arg_parser():
    import argparse
    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program merge a series of ProtoMS template files")
    parser.add_argument(
        '-f', '--files', nargs="+", help="the name of the template files")
    parser.add_argument(
        '-o', '--out', help="the name of the merged template file")
    return parser


if __name__ == "__main__":

    args = get_arg_parser().parse_args()

    # Setup the logger
    logger = simulationobjects.setup_logger("merge_templates_py.log")

    if args.files is None:
        raise simulationobjects.SetupError(
            "Nothing to do! No template files provided. Exiting.")
    if args.out is None:
        raise simulationobjects.SetupError(
            "A name for the merged template file needs to be provided!")

    tem = merge_templates(args.files)
    if isinstance(tem, simulationobjects.TemplateFile):
        tem.write(args.out)
