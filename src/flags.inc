
c   This file holds the flags the control the simulation

c     whether or not we have a solvent box or cap
      logical SolventBox,SolventCap
c     whether we will use the parameters in the protein pdb
      logical PdbParams
c     whether we are doing a jaws1 or jaws2 simulation
      logical doJaws
c     whether the jaws2 results for each move will be printed to results
      logical everyjaws2
c     whether we will print the free energy data in results
      logical printFE, mbar
c     whether or not to neutralise the protein
      logical NeutraliseProtein
c     whether or not to print debugging output
      logical protodebug
c     whether or not we are doing full energy calculations (for testing)
      logical testEnergy
c     whether or not an initial TotalEnergy has been calculated
      logical GotInitialTotalEnergy
c     whether or not lambdaf/lambdab are different to lambda
      logical DifferentLambdas
c     whether we will do single-topology with lambda free energy calculation
      logical SingleTopologyFreeNrg
c     whether we will do dual-topology with lambda free energy calculation
      logical DualTopologyFreeNrg
c     whether we are using Lorentz-Berthelot combining rules
      logical LorentzBerthelot
c     whether or not to do preferential sampling
      logical PrefSample
c     whether or not to do LambdaRE and TemperatureRE
      logical doLambdaRE,doTemperatureRE,doREST,doMultiGCMC,doMultiJaws2
c     whether or not to calculate lambda-gradient and not energy
      logical doDerivative,doLowDerivative
c     whether or not the MBAR, BAR and TI data needs to be re-calculated
      logical reCalculateFreeEnergies
c     whether or not to perform a dihedral flexibilty tuning simulation
      logical tuneDihedral

      common /LSIMFLAGS/ SolventBox,SolventCap,NeutraliseProtein,
     .          protodebug,
     .          GotInitialTotalEnergy,testEnergy,
     .          SingleTopologyFreeNrg,DualTopologyFreeNrg,
     .          DifferentLambdas,LorentzBerthelot,PrefSample,
     .          doLambdaRE,doTemperatureRE,doREST,doMultiGCMC,doMultiJaws2,
     .          PdbParams,doJaws,everyjaws2,
     .          printFE,mbar,doDerivative,doLowDerivative,reCalculateFreeEnergies,
     .		tuneDihedral

      save /LSIMFLAGS/

      integer cuttype
      integer boundary
      common /ISIMFLAGS/ cuttype,boundary
      
      save /ISIMFLAGS/
