
c   This file holds the variables that store the simulation parameters (temperature etc.)

c     the current value of lambda, and the forwards and backwards values (FEP)
      double precision Lambda,LambdaB,LambdaF
c     the temperature in C and K
      double precision temperatureC,temperature 
c     the charge you want on the protein
      double precision desiredCharge
c     B factor
      double precision B     
c     bias-jaws
      double precision bias  
c     one over the default coulombic 1-4 scale factor
      double precision OneOverScl14C
c     one over the default LJ 1-4 scale factor
      double precision OneOverScl14LJ
c     parameters associated with spherical cap
      double precision cap(3),capradius,oldcapradius,capconstant
c     space for the coordinates of 4 atoms (used for evaluating internal
c     angles of atoms)
      double precision TempAtoms(4,3)
c     solvent box dimensions...
      double precision BoxDimension(3),BoxDimension2(3)
      double precision BoxOrigin(3),BoxTop(3)
      double precision OldBoxDimension(3),OldBoxOrigin(3)

  
c     grid origin for GCMC
      double precision GCtheta,jaws1,jaws2,thres
      double precision coordx
      double precision coordy
      double precision coordz
c     grid dimensions for GCMC
      double precision gridx
      double precision gridy
      double precision gridz
      double precision space
c     grid center for GCMC/JAWS
      double precision centerx
      double precision centery
      double precision centerz
c     grid origin for GCMC
      double precision coordhwx
      double precision coordhwy
      double precision coordhwz
c     grid dimensions for GCMC
      double precision gridhwx
      double precision gridhwy
      double precision gridhwz


c     the size of the non-bonded cutoff (and squared version)
      double precision cut,cut2
c     the size of the feathering of the non-bonded cutoff (and squared)
      double precision feather,feather2,featherfac
c     JM The parameters of the softcore
      integer CoulSoftPower,LJSoftPower,GBSoftPower
      double precision DeltaCore
c     SG Additional softcore parameters
      double precision DeltaCoreCoul,DeltaCoreCoulL,SoftPower,SoftPower2,SoftPowerInv6,SoftPowerInv

c     different soft-core parameters
      integer SOFTOLD,SOFT66,SOFTAMB,SoftType
      parameter (SOFTOLD=1)
      parameter (SOFT66=2)
      parameter (SOFTAMB=3)

c     solute-tempering
      double precision resttempC,resttemp,rest_beta
      double precision restfac1,restfac2
c     restscaling(1) : intra-solute bond-angle
c     restscaling(2) : intra-solute dihedral
c     restscaling(3) : intra-solute lj
c     restscaling(4) : intra-solute coul
c     restscaling(5) : solute-solute
c     restscaling(6) : solute-protein
c     restscaling(7) : solute-water
      double precision restscaling(7)
      integer restlevel(7)
      
c     whether or not to calculate energies in the same way as Sire
c     (this slows things down, but is useful to allow debugging of Sire)
      logical SIRE_COMPATIBILITY_MODE 

      common /SIMPARAMVALS/ cut,cut2,centerx,centery,centerz,
     .           Lambda,LambdaF,LambdaB,BoxDimension,BoxOrigin,feather,feather2,featherfac,
     .           OldBoxDimension,OldBoxOrigin,temperatureC,temperature,B,coordx,coordy,coordz,gridx,gridy,gridz,space,bias,
     .           BoxDimension2,BoxTop,TempAtoms,desiredCharge,coordhwx,coordhwy,coordhwz,gridhwx,gridhwy,gridhwz,
     .           cap,capradius,oldcapradius,capconstant,GCtheta,jaws1,jaws2,
     .           OneOverScl14C,OneOverScl14LJ,thres,
     .             SIRE_COMPATIBILITY_MODE  
      common /SIMPARAMVALS2/ DeltaCoreCoul,DeltaCoreCoulL,SoftPower,SoftPowerInv,SoftPower2,SoftPowerInv6,
     .                          DeltaCore,CoulSoftPower,LJSoftPower,GBSoftPower
      
      save /SIMPARAMVALS/
      save /SIMPARAMVALS2/
      
      common /REST/ restscaling,restfac1,restfac2,resttempC,resttemp,rest_beta,restlevel
      save /REST/
