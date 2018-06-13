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

from . import simulationobjects

logger = logging.getLogger('protoms')


def merge_templates(templates):
    """
    Merge a series of ProtoMS template files

    Parameters
    ----------
    templates : list of string
      the names of the template file

    Returns
    -------
    TemplateFile
      the merged template file
    or
      string
    the filename of the single unique template
    """

    logger.debug("Running merge_templates with arguments: ")
    logger.debug("\ttemplates = %s" % " ".join(templates))
    logger.debug(
        "This will merge all templates, renumbering force field parameters")

    # Make it a unique list
    templates2 = []
    for t in templates:
        if t not in templates2:
            templates2.append(t)

    if len(templates2) == 1:
        return templates2[0]

    temfile = simulationobjects.TemplateFile(templates2[0])
    for t in templates2[1:]:
        temfile2 = simulationobjects.TemplateFile(t)
        temfile.append(temfile2)
    return temfile


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
