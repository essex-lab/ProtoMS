# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

"""
Program to produce a template file for a solute in the format
required to proceed with a jaws1 simulation for the fragments.

This new template can be based on a previously given standard
template for a solute, or produced from scratch from a pdb file
of a solute.

Part of the parameters required to be added to the template file will
be taken from the output folder of a hydration free energy calculation,
required to have been run before the jaws for the fragments.

Some of the parameters however will be user provided (though some
defaults are stablished).

Can be executed from the command line as a stand-alone program
"""

import simulationobjects

import calc_ti

def modify_template_jfrag(tempin,jtheta=[0.15],jcorr=[1.0],jpmf=[]) :
  """
  Writes out a template for simulations of jaws1 for the fragments
  from a previous standard solute template

  Parameters
  ----------
  tempin : string
    the file name of the original solute template
  jtheta : list of floats
    the range of a theta move for the each solute
  jcorr : list of float
    the concentration correction parameter to include
    in the template for each of the solutes
  jpmf : list of lists or list of numpy arrays
    the indeces of the polynomial function that best fits
    the hydration free energy calculation for each solute

  Returns
  -------
  template object
    the template object where the new template is stored
    printed to
  """

  tempfile = simulationobjects.TemplateFile(filename=tempin)
  for ind,moltemplate in enumerate(tempfile.templates) :
    moltemplate.jtheta = jtheta[ind]
    moltemplate.jcorr = jcorr[ind]
    moltemplate.jpmf = jpmf[ind]
  return tempfile

if __name__ == '__main__' :

  import argparse

  parser = argparse.ArgumentParser(description="Produce a jaws-fragment template. You need to have run previously the hydration free energy simulation of this molecule.")
  parser.add_argument('-t','--template',help="The standard solute template (.tem) file of your molecule. It was necesarelly produced for your hydration free energy calculation.")
  parser.add_argument('-o','--outtemplate',help="The name of the file where the jaws-fragment template will be writen to. Default = fragment.tem",default='fragment.tem')
  parser.add_argument('-f','--hydfolder',help="The root directory that contains all the output files of the hydration free energy simulation. Default = out",default="out")
  parser.add_argument('-r','--results',help="The name of the results file in the hydration free energy calculation. Default = results.",default='results')
  parser.add_argument('-jt','--jtheta',help="The maximum amplitud of the solute-theta moves. Default = 0.15",default=0.15)
  parser.add_argument('-jc','--jcorr',help="The concentration correction to be aplied in the jaws simulation. Default = 1.0",default= 1.0)
  parser.add_argument('--skip',help="The number of snapshots to skip in the calculation of the hydration free energy. Default = 0",default=0)
  parser.add_argument('--max',help="The highest snapshot to use in the calculation of the hdyration free energy. The default should use all. Deafult = 99999",default=99999)
  parser.add_argument('--minorder',help="The minimum order of the polynomial to fit the PMF of the hydration free energy. Default = 4",default=4)
  parser.add_argument('--maxorder',help="The maximum order of the polynomial to fit the PMF of the hydration free energy. Default = 5",default=5)
  parser.add_argument('--fitname',help="The name of the file where the fit of the hydration PMF will be saved. Deafult=fit.png",default='fit.png')
  args = parser.parse_args()

  # Setup the logger
  logger = simulationobjects.setup_logger()

  if args.hydfolder[-1] != "/" : args.hydfolder = args.hydfolder + "/"

  verbose = {"total":False,"gradient":False,"pmf":False,"uncert":False,"lambda":False}  
  
  num_grad_kind = 'both'
  anal_grad = True
  fit = True

  lambdas,gradients,grad_std,pmf,pmf_std = calc_ti.ti(args.hydfolder,args.results,args.skip,args.max,verbose,num_grad_kind,anal_grad)

  fit_coef = calc_ti.fit_pmf(lambdas,pmf,int(args.minorder),int(args.maxorder),args.fitname)

  newtemplate = modify_template_jfrag(args.template,jtheta=[float(args.jtheta)],jcorr=[float(args.jcorr)],jpmf=[fit_coef])
  newtemplate.write(filename=args.outtemplate)
  logger.info("The jaws-fragment template file has been produced in %s.\n"%args.outtemplate)

  


























    
