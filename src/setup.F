      subroutine setup
c###########################################################
c
c     This routine sets up the system, loads molecules
c     etc.
c
c     (C) Christopher Woods      9 August 2002
c
c###########################################################

c     call a routine to get all the parameters to the program 
c     - this loads all of the parameter files and pdb files
      call getOptions

c     now the system is loaded, cache all of the variable groups
      call scanForVariableGroups

c     now precalculate all the pairs of charge/LJ (clj) parameters
c     used in this simulation
      call calculateCLJPairs

      end
