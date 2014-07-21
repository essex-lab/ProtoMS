# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross


""" 
Routines that wrap AmberTools programs

This module defines two public functions:
run_antechamber
run_parmchk 
    
Can be executed from the command line as a stand-alone program
"""

import os
import subprocess
import tempfile

import simulationobjects

def _run_program ( name, command ):
    """ 
    Wrapper for an AmberTools program with default parameters

    run_program ( name, command )

    Parameters
    ----------        
    name : string
    	the name of the program to run
    command : string
    	the command to execute it
    	
    Raises
    ------
    SetupError
    	if the program failed to execute properly
    """

    # Create temporary file to write stdout, stderr of the program
    tmpfile,tmpname = tempfile.mkstemp()
    ret_code = subprocess.call ( command, shell = True , stdout=tmpfile, stderr=tmpfile)
    # Catch some error codes
    if ret_code == 127:
        raise simulationobjects.SetupError ( "Unable to find %s executable, please make sure this is present in your PATH."%name )
    if ret_code == 1:
        # Get the error message from the temporary file
        errmsg = "\n".join(line for line in open(tmpname).readlines())
        os.remove(tmpname)
        raise simulationobjects.SetupError ( "%s was not able to run successfully. Please check output. Error message was:\n%s"%(name,errmsg)  )

    os.remove(tmpname)

def run_antechamber ( lig, charge, resnam = None):
    """ 
    Wrapper for antechamber with default parameters and AM1-BCC charges

    Parameters
    ----------        
    lig : PDBFile or string
    	the ligand to run Antechamber on 
    charge : int
    	the net charge of the ligand
    resnam : string, optional
    	the residue name of the ligand
    	
    Returns
    -------
    string
    	the filename of the created prepi file
    """

    # Check if it is a string, otherwise assumes it is a pdb object
    if isinstance(lig,basestring) :
      name = lig
    else :
      name = lig.name

    if resnam is None :
      resnamstr = ""
    else :
      resnamstr = "-rn "+resnam

    # Remove the extension from the filename
    out_name = os.path.splitext ( name )[0]
    cmd = 'antechamber -i %s -fi pdb -o %s.prepi -fo prepi -c bcc -nc %d %s -pf y' % ( name, out_name,  charge, resnamstr )
    _run_program ( "antechamber", cmd)
    subprocess.call ( "rm ANTECHAMBER_* ATOMTYPE.INF PREP.INF NEWPDB.PDB sqm.in sqm.out sqm.pdb", shell = True )
    return "%s.prepi"%out_name

def run_parmchk ( lig ):
    """ 
    Wrapper for parmcheck with default parameters

    Parameters
    ----------        
    lig : PDBFile or string
		the ligand to run Parmcheck on
		
    Returns
    -------
    string
    	the filename of the created frcmod file
    """
   
    # Check if it is a string, otherwise assumes it is a pdb object
    if isinstance(lig,basestring) :
      name = lig
    else :
      name = lig.name

    # Remove the extension from the filename
    out_name = os.path.splitext ( name )[0]
    cmd = 'parmchk -i %s.prepi -f prepi -o %s.frcmod' % ( out_name, out_name )
    _run_program ( "parmchk", cmd)
    return "%s.frcmod"%out_name

if __name__ == "__main__":

    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(description="Program to run antechamber and parmchk for a series of PDB-files")
    parser.add_argument('-f','--files',nargs="+",help="the name of the PDB-files")
    parser.add_argument('-n','--name',help="the name of the solute",default="UNK")
    parser.add_argument('-c','--charge',nargs="+",type=float,help="the net charge of each PDB-file")
    args = parser.parse_args()

    for i,filename in enumerate(args.files) :
      if args.charge is None or i >= len(args.charge) :
        charge = 0.0
      else :
        charge = args.charge[i]
      run_antechamber(filename,charge,args.name)
      run_parmchk(filename+"P")
