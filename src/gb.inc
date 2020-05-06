c
c     File created by Julien Michel 25/06/04
c     Include files for Generalised Born calculations with ProtoMC
c

c     Whether or not to use the GB model
      logical doGB
c     Maximum number of GB types
      integer MAXGBTYPES
      parameter(MAXGBTYPES=200)
c     Conversion from SIGMA to RSTAR
      double precision SIGTORS
      parameter(SIGTORS=0.5612310242D+00)
c     The minimum size of an intrinsic born radius
      double precision MINIBORN
      parameter(MINIBORN=1.100+00)
      integer PDA
      parameter (PDA=1)
c     These Arrays store GB parameters
      character*2 GBType(MAXGBTYPES)
      logical GotGB(MAXGBTYPES,2)
      double precision GBIBorn(MAXGBTYPES)
      double precision GBScale(MAXGBTYPES)
      double precision IntrBornDB(MAXCLJ)
      double precision PDAScaleDB(MAXCLJ)
c    This is for perturbed atoms
      double precision VariableIntrBorn(MAXVARIABLEPARAM,5)
      double precision VariablePDAScale(MAXVARIABLEPARAM,5)
c
c     2D Arrays that holds the Born radii and scale factors of the PDA     
c     for Solutes, Proteins and Solvents
c      double precision ZSolIBorn(MAXSOLUTES,MAXSOLUTEATOMS) 
c      double precision ZSolScaleBorn(MAXSOLUTES,MAXSOLUTEATOMS) 
      double precision ZSolContBorn(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision NewZSolContBorn(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision OldZSolContBorn(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision StackZSolContBorn(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision ZSolEffBorn(MAXSOLUTES,MAXSOLUTEATOMS) 
c      double precision ZProtIBorn(MAXPROTEINS,MAXATOMS)  
c      double precision ZProtScaleBorn(MAXPROTEINS,MAXATOMS) 
      double precision ZProtContBorn(MAXPROTEINS,MAXATOMS)
      double precision NewZProtContBorn(MAXPROTEINS,MAXATOMS)
      double precision OldZProtContBorn(MAXPROTEINS,MAXATOMS)
      double precision ZProtContBornF(MAXPROTEINS,MAXATOMS)
      double precision NewZProtContBornF(MAXPROTEINS,MAXATOMS)
      double precision OldZProtContBornF(MAXPROTEINS,MAXATOMS)
      double precision ZProtContBornB(MAXPROTEINS,MAXATOMS)
      double precision NewZProtContBornB(MAXPROTEINS,MAXATOMS)
      double precision OldZProtContBornB(MAXPROTEINS,MAXATOMS)
      double precision StackZProtContBorn(MAXPROTEINS,MAXATOMS)
      double precision StackZProtContBornF(MAXPROTEINS,MAXATOMS)
      double precision StackZProtContBornB(MAXPROTEINS,MAXATOMS)
      double precision ZProtEffBorn(MAXPROTEINS,MAXATOMS) 
      double precision ZProtEffBornF(MAXPROTEINS,MAXATOMS) 
      double precision ZProtEffBornB(MAXPROTEINS,MAXATOMS) 
c      double precision ZSvnIBorn(MAXSOLVENTS,MAXSOLVENTATOMS) 
c      double precision ZSvnScaleBorn(MAXSOLVENTS,MAXSOLVENTATOMS) 
      double precision StackZSvnContBorn(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision StackZSvnContBornF(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision StackZSvnContBornB(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision ZSvnContBorn(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision NewZSvnContBorn(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision OldZSvnContBorn(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision ZSvnContBornF(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision NewZSvnContBornF(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision OldZSvnContBornF(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision ZSvnContBornB(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision NewZSvnContBornB(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision OldZSvnContBornB(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision ZSvnEffBorn(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision ZSvnEffBornF(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision ZSvnEffBornB(MAXSOLVENTS,MAXSOLVENTATOMS)
c Cutoff for the computation of the effective Born radius
      double precision BornCut
      double precision BornCut2
c    For FEP
      double precision ZSolEffBornF(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision ZSolEffBornB(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision ZSolContBornF(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision NewZSolContBornF(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision OldZSolContBornF(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision ZSolContBornB(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision NewZSolContBornB(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision OldZSolContBornB(MAXSOLUTES,MAXSOLUTEATOMS)            
      double precision StackZSolContBornF(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision StackZSolContBornB(MAXSOLUTES,MAXSOLUTEATOMS)

      common /DGBORN/ ZSolEffBorn,ZSolContBorn,ZProtEffBorn,
     .    ZProtContBorn,ZSvnEffBorn,ZSvnContBorn,
     .    GBIBorn,GBScale,IntrBornDB,PDAScaleDB,
     .  NewZSolContBorn,NewZProtContBorn,NewZSvnContBorn,
     .  OldZSolContBorn,OldZProtContBorn,OldZSvnContBorn,
     .    ZProtContBornF,NewZProtContBornF,OldZProtContBornF,
     .    ZProtContBornB,NewZProtContBornB,OldZProtContBornB,
     .    StackZSolContBorn,StackZProtContBorn,StackZSvnContBorn,
     .    BornCut,BornCut2,VariableIntrBorn,VariablePDAScale,
     .    ZSolEffBornF,ZSolEffBornB,NewZSolContBornF,OldZSolContBornF,
     .    NewZSolContBornB,OldZSolContBornB,ZSolContBornF,ZSolContBornB,
     .    StackZSolContBornF,StackZSolContBornB,
     .    ZSvnContBornF,NewZSvnContBornF,OldZSvnContBornF,
     .    ZSvnContBornB,NewZSvnContBornB,OldZSvnContBornB,ZProtEffBornF,
     .    ZProtEffBornB,ZSvnEffBornF,ZSvnEffBornB,StackZProtContBornF,
     .  StackZProtContBornB,StackZSvnContBornF,StackZSvnContBornB

      common /LGBORN/ doGB,GotGB

      common /CGBORN/ GBType

      save /LGBORN/
      save /CGBORN/
c
c Arrays that holds the GB energy components
c 
      double precision OldGBEnergy
      double precision OldGBEnergyF
      double precision OldGBEnergyB
      double precision NewGBEnergy
      double precision NewGBEnergyF
      double precision NewGBEnergyB
      double precision AvgGBEnergy
      double precision AvgGBEnergyF
      double precision AvgGBEnergyB
      double precision AvgGBEnergy2
      double precision AvgGBEnergyF2
      double precision AvgGBEnergyB2

      common /DGBNRG/ OldGBEnergy,OldGBEnergyF,OldGBEnergyB,
     .  NewGBEnergy,NewGBEnergyF,NewGBEnergyB,
     .    AvgGBEnergy,AvgGBEnergyF,AvgGBEnergyB,
     .    AvgGBEnergy2,AvgGBEnergyF2,AvgGBEnergyB2

      save /DGBNRG/

c     Whether or not the GB energy should be recomputed for one given atom
      logical SolAtomdoBorn(MAXSOLUTES,MAXSOLUTEATOMS)
      logical ProtAtomdoBorn(MAXPROTEINS,MAXATOMS)
      logical SvnAtomdoBorn(MAXSOLVENTS,MAXSOLVENTATOMS)
      logical calcResidueGB(MAXPROTEINS,MAXRESIDUES)
      logical calcSoluteGB(MAXSOLUTES)
      logical calcSolventGB(MAXSOLVENTS)
c     If this residue has just been moved
      integer MovedResidue(MAXPROTEINS,MAXRESIDUES)
c     If this solute has just been moved
      integer MovedSolute(MAXSOLUTES)
c     If this solvent has just been moved
      integer MovedSolvent(MAXSOLVENTS)
c     Recompute born radii if they changed by more than gbthreshold
      double precision gbthreshold
c     Store old born radii
      double precision OldZSolEffBorn(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision OldZProtEffBorn(MAXPROTEINS,MAXATOMS)
      double precision OldZSvnEffBorn(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision OldZSolEffBornF(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision OldZProtEffBornF(MAXPROTEINS,MAXATOMS)
      double precision OldZSvnEffBornF(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision OldZSolEffBornB(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision OldZProtEffBornB(MAXPROTEINS,MAXATOMS)
      double precision OldZSvnEffBornB(MAXSOLVENTS,MAXSOLVENTATOMS)

c     To be able to save coordinates of moving things...
      integer countGBSolute
      integer countGBResidue
      integer countGBSolvent
      integer gbSoluteList(MAXSOLUTES)
      integer gbResidueList(MAXRESIDUES)
      integer gbResProtList(MAXRESIDUES)
      integer gbSolventList(MAXSOLVENTS)
      double precision gbSoluteCoords(MAXSOLUTES,MAXSOLUTEATOMS,3)
      double precision gbSoluteCoordsF(MAXSOLUTES,MAXSOLUTEATOMS,3)
      double precision gbSoluteCoordsB(MAXSOLUTES,MAXSOLUTEATOMS,3)
      double precision gbSolCentre(MAXSOLUTES,3)
      double precision gbSolRadius(MAXSOLUTES)
      double precision gbSolventCoords(MAXSOLVENTS,MAXSOLVENTATOMS,3)
      double precision gbSvnCentre(MAXSOLVENTS,3)
      integer gbResidueIndex(MAXPROTEINS,MAXATOMS)
      double precision gbResidueCoords(MAXRESIDUES,MAXSCATOMS,3)
      double precision gbResCentre(MAXRESIDUES,3)
      double precision gbResRadius(MAXRESIDUES)

      common /DFASTGB/ OldZSolEffBorn,OldZProtEffBorn,OldZSvnEffBorn,
     .    OldZSolEffBornF,OldZProtEffBornF,OldZSvnEffBornF,
     .    OldZSolEffBornB,OldZProtEffBornB,OldZSvnEffBornB,
     .  gbthreshold,
     .    gbSoluteCoords,gbSoluteCoordsF,gbSoluteCoordsB,
     .    gbSolCentre,gbSolRadius,
     .    gbResidueCoords,gbResCentre,gbResRadius,
     .    gbSolventCoords,gbSvnCentre

      common /LFASTGB/ SolAtomdoBorn,ProtAtomdoBorn,SvnAtomdoBorn,
     . calcResidueGB,calcSoluteGB,calcSolventGB

      common /IFASTGB/ MovedResidue,MovedSolute,MovedSolvent,
     .    countGBSolute,countGBResidue,countGBSolvent,
     .    gbSoluteList,gbResidueList,gbSolventList,gbResProtList,
     .  gbResidueIndex

      save /DFASTGB/
      save /LFASTGB/
      save /IFASTGB/

c    Stuff for the smoothing function of Onufriev Bashford Case
c    this is GB OBC(II) here
      double precision ONUA
      parameter(ONUA=1.0d0)
      double precision ONUB
      parameter(ONUB=-0.80d0) 
      double precision ONUG
      parameter(ONUG=4.85d0)
      logical gbproteins

      common /LONUFGB/ gbproteins
      save /LONUFGB/

