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

import simulationobjects

logger = logging.getLogger('protoms')

def merge_templates(templates) :
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
  """
  
  logger.debug("Running merge_templates with arguments: ")
  logger.debug("\ttemplates = %s"%" ".join(templates)) 
  logger.debug("This will merge all templates, renumbering force field parameters")
  
  templates = list(set(templates)) # Make it a unique list
  temfile = simulationobjects.TemplateFile(templates[0])
  for t in templates[1:] :
    temfile2 = simulationobjects.TemplateFile(t)
    temfile.append(temfile2)
  return temfile

if __name__ == "__main__":

  import argparse
  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program merge a series of ProtoMS template files")
  parser.add_argument('-f','--files',nargs="+",help="the name of the template files")
  parser.add_argument('-o','--out',help="the name of the merged template file")
  args = parser.parse_args()

  # Setup the logger
  logger = simulationobjects.setup_logger("merge_templates_py.log")

  if args.files is None : 
    print "Nothing to do! Exiting."
  
  tem = merge_templates(args.files)
  tem.write(args.out)
