c   This file holds the stuff that controls the use of a quick sampling
c   potential
c     Energies in state i and j of accurate potential
      integer ACCNRGI,ACCNRGJ,QUICKNRGI,QUICKNRGJ
      parameter(ACCNRGI=1)
      parameter(ACCNRGJ=2)
      parameter(QUICKNRGI=3)
      parameter(QUICKNRGJ=4)
      integer STOREDNRGS
      parameter(STOREDNRGS=4)
      integer CONFI,CONFJ
      parameter(CONFI=1)
      parameter(CONFJ=2)
      integer STOREDCONFS
      parameter(STOREDCONFS=2)
c
c     Variables storing the energy of one configuration
c
      double precision StoredProBnd(MAXPROTEINS,STOREDNRGS)
      double precision StoredProAng(MAXPROTEINS,STOREDNRGS)
      double precision StoredProDih(MAXPROTEINS,STOREDNRGS)
      double precision StoredProUBy(MAXPROTEINS,STOREDNRGS)
      double precision StoredProNBC(MAXPROTEINS,STOREDNRGS)
      double precision StoredProNBLJ(MAXPROTEINS,STOREDNRGS)
      double precision StoredProBndF(MAXPROTEINS,STOREDNRGS)
      double precision StoredProAngF(MAXPROTEINS,STOREDNRGS)
      double precision StoredProDihF(MAXPROTEINS,STOREDNRGS)
      double precision StoredProUByF(MAXPROTEINS,STOREDNRGS)
      double precision StoredProNBCF(MAXPROTEINS,STOREDNRGS)
      double precision StoredProNBLJF(MAXPROTEINS,STOREDNRGS)
      double precision StoredProBndB(MAXPROTEINS,STOREDNRGS)
      double precision StoredProAngB(MAXPROTEINS,STOREDNRGS)
      double precision StoredProDihB(MAXPROTEINS,STOREDNRGS)
      double precision StoredProUByB(MAXPROTEINS,STOREDNRGS)
      double precision StoredProNBCB(MAXPROTEINS,STOREDNRGS)
      double precision StoredProNBLJB(MAXPROTEINS,STOREDNRGS)
      double precision StoredProProC(MAXPROTEINS,MAXPROTEINS,STOREDNRGS)
      double precision StoredProProLJ(MAXPROTEINS,MAXPROTEINS,STOREDNRGS)
      double precision StoredProProCF(MAXPROTEINS,MAXPROTEINS,STOREDNRGS)
      double precision StoredProProLJF(MAXPROTEINS,MAXPROTEINS,STOREDNRGS)
      double precision StoredProProCB(MAXPROTEINS,MAXPROTEINS,STOREDNRGS)
      double precision StoredProProLJB(MAXPROTEINS,MAXPROTEINS,STOREDNRGS)
      double precision StoredProProDist(MAXPROTEINS,MAXPROTEINS,STOREDNRGS)

      double precision StoredSolBnd(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolAng(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolDih(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolUBy(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolNBC(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolNBLJ(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolBndF(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolAngF(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolDihF(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolUByF(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolNBCF(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolNBLJF(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolBndB(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolAngB(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolDihB(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolUByB(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolNBCB(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolNBLJB(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolSolC(MAXSOLUTES,MAXSOLUTES,STOREDNRGS)
      double precision StoredSolSolLJ(MAXSOLUTES,MAXSOLUTES,STOREDNRGS)
      double precision StoredSolSolCF(MAXSOLUTES,MAXSOLUTES,STOREDNRGS)
      double precision StoredSolSolLJF(MAXSOLUTES,MAXSOLUTES,STOREDNRGS)
      double precision StoredSolSolCB(MAXSOLUTES,MAXSOLUTES,STOREDNRGS)
      double precision StoredSolSolLJB(MAXSOLUTES,MAXSOLUTES,STOREDNRGS)
      double precision StoredSolSolDist(MAXSOLUTES,MAXSOLUTES,STOREDNRGS)
      double precision StoredProSolC(MAXPROTEINS,MAXSOLUTES,STOREDNRGS)
      double precision StoredProSolLJ(MAXPROTEINS,MAXSOLUTES,STOREDNRGS)
      double precision StoredProSolCF(MAXPROTEINS,MAXSOLUTES,STOREDNRGS)
      double precision StoredProSolLJF(MAXPROTEINS,MAXSOLUTES,STOREDNRGS)
      double precision StoredProSolCB(MAXPROTEINS,MAXSOLUTES,STOREDNRGS)
      double precision StoredProSolLJB(MAXPROTEINS,MAXSOLUTES,STOREDNRGS)
      double precision StoredProSolDist(MAXPROTEINS,MAXSOLUTES,STOREDNRGS)

      double precision StoredSvnSvnC(STOREDNRGS)
      double precision StoredSvnSvnLJ(STOREDNRGS)
      double precision StoredSvnSvnCF(STOREDNRGS)
      double precision StoredSvnSvnLJF(STOREDNRGS)
      double precision StoredSvnSvnCB(STOREDNRGS)
      double precision StoredSvnSvnLJB(STOREDNRGS)
      double precision StoredSvnSvnDist(STOREDNRGS)
      double precision StoredProSvnC(MAXPROTEINS,STOREDNRGS)
      double precision StoredProSvnLJ(MAXPROTEINS,STOREDNRGS)
      double precision StoredProSvnCF(MAXPROTEINS,STOREDNRGS)
      double precision StoredProSvnLJF(MAXPROTEINS,STOREDNRGS)
      double precision StoredProSvnCB(MAXPROTEINS,STOREDNRGS)
      double precision StoredProSvnLJB(MAXPROTEINS,STOREDNRGS)
      double precision StoredProSvnDist(MAXPROTEINS,STOREDNRGS)
      double precision StoredSolSvnC(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolSvnLJ(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolSvnCF(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolSvnLJF(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolSvnCB(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolSvnLJB(MAXSOLUTES,STOREDNRGS)
      double precision StoredSolSvnDist(MAXSOLUTES,STOREDNRGS)
      double precision StoredSvnProb(MAXSOLVENTS,STOREDNRGS)
      double precision StoredSumProb(STOREDNRGS)
      double precision StoredMaxProb(STOREDNRGS)
c
c     The extra energies...
c
      double precision StoredSolventCap(STOREDNRGS)
      double precision StoredSolventCapPart(STOREDNRGS)
      double precision StoredHarmEnergy(STOREDNRGS)
      double precision StoredHarmEnergyF(STOREDNRGS)
      double precision StoredHarmEnergyB(STOREDNRGS)
      double precision StoredHarmEnergyPart(STOREDNRGS)
      double precision StoredHarmEnergyPartF(STOREDNRGS)
      double precision StoredHarmEnergyPartB(STOREDNRGS)
      double precision StoredFBotEnergy(STOREDNRGS)
      double precision StoredFBotEnergyF(STOREDNRGS)
      double precision StoredFBotEnergyB(STOREDNRGS)
      double precision StoredFBotEnergyPart(STOREDNRGS)
      double precision StoredFBotEnergyPartF(STOREDNRGS)
      double precision StoredFBotEnergyPartB(STOREDNRGS)
      double precision StoredTotalSurface(STOREDNRGS)
      double precision StoredTotalSurfaceF(STOREDNRGS)
      double precision StoredTotalSurfaceB(STOREDNRGS)
      double precision StoredSAEnergy(STOREDNRGS)
      double precision StoredSAEnergyF(STOREDNRGS)
      double precision StoredSAEnergyB(STOREDNRGS)
      double precision StoredGBEnergy(STOREDNRGS)
      double precision StoredGBEnergyF(STOREDNRGS)
      double precision StoredGBEnergyB(STOREDNRGS)
      double precision StoredExtraEnergyPart(STOREDNRGS)
      double precision StoredExtraEnergyPartF(STOREDNRGS)
      double precision StoredExtraEnergyPartB(STOREDNRGS)
      double precision StoredExtraEnergy(STOREDNRGS)
      double precision StoredExtraEnergyF(STOREDNRGS)
      double precision StoredExtraEnergyB(STOREDNRGS)

      double precision StoredTotalEnergy(STOREDNRGS)
      double precision StoredTotalEnergyF(STOREDNRGS)
      double precision StoredTotalEnergyB(STOREDNRGS)



      common /DSTOREATA/ StoredProBnd,StoredProAng,StoredProDih,StoredProUBy,
     . StoredProNBC,StoredProNBLJ,StoredProBndF,StoredProAngF,StoredProDihF,
     . StoredProUByF,StoredProNBCF,StoredProNBLJF,StoredProBndB,StoredProAngB,
     . StoredProDihB,StoredProUByB,StoredProNBCB,StoredProNBLJB,StoredProProC,
     . StoredProProLJ,StoredProProCF,StoredProProLJF,StoredProProCB,StoredProProLJB,
     . StoredProProDist,StoredSolBnd,StoredSolAng,StoredSolDih,StoredSolUBy,
     . StoredSolNBC,StoredSolNBLJ,StoredSolBndF,StoredSolAngF,StoredSolDihF,
     . StoredSolUByF,StoredSolNBCF,StoredSolNBLJF,StoredSolBndB,StoredSolAngB,
     . StoredSolDihB,StoredSolUByB,StoredSolNBCB,StoredSolNBLJB,StoredSolSolC,
     . StoredSolSolLJ,StoredSolSolCF,StoredSolSolLJF,StoredSolSolCB,StoredSolSolLJB,
     . StoredSolSolDist,StoredProSolC,StoredProSolLJ,StoredProSolCF,StoredProSolLJF,
     . StoredProSolCB,StoredProSolLJB,StoredProSolDist,StoredSvnSvnC,StoredSvnSvnLJ,
     . StoredSvnSvnCF,StoredSvnSvnLJF,StoredSvnSvnCB,StoredSvnSvnLJB,StoredSvnSvnDist,
     . StoredProSvnC,StoredProSvnLJ,StoredProSvnCF,StoredProSvnLJF,StoredProSvnCB,
     . StoredProSvnLJB,StoredProSvnDist,StoredSolSvnC,StoredSolSvnLJ,StoredSolSvnCF,
     . StoredSolSvnLJF,StoredSolSvnCB,StoredSolSvnLJB,StoredSolSvnDist

      common /DSTORETWODATA/ StoredSvnProb,StoredSumProb,StoredMaxProb,
     . StoredTotalEnergy,StoredTotalEnergyF,StoredTotalEnergyB,StoredSolventCap,
     . StoredHarmEnergy,StoredHarmEnergyF,StoredHarmEnergyB,StoredFBotEnergy,
     . StoredFBotEnergyF,StoredFBotEnergyB,StoredTotalSurface,StoredTotalSurfaceF,
     . StoredTotalSurfaceB,StoredSAEnergy,StoredSAEnergyF,StoredSAEnergyB,StoredGBEnergy,
     . StoredGBEnergyF,StoredGBEnergyB,StoredExtraEnergyPart,StoredExtraEnergyPartF,
     . StoredExtraEnergyPartB,StoredExtraEnergy,StoredExtraEnergyF,
     . StoredExtraEnergyB,StoredFBotEnergyPart,StoredFBotEnergyPartF,
     . StoredFBotEnergyPartB,StoredHarmEnergyPart,StoredHarmEnergyPartF,
     . StoredHarmEnergyPartB,StoredSolventCapPart

c    Keeping track of saved coordinates
      double precision StoreZCoords(MAXPROTEINS,MAXATOMS,3,STOREDCONFS)
      double precision StoreZSolCoords(MAXSOLUTES,MAXSOLUTEATOMS,3,STOREDCONFS)
      double precision StoreSolCentre(MAXSOLUTES,3,STOREDCONFS)
      double precision StoreSolGeomCentre(MAXSOLUTES,3,STOREDCONFS)
      double precision StoreSolDummies(MAXSOLUTES,3,3,STOREDCONFS)
      double precision StoreZSvnCoords(MAXSOLVENTS,MAXSOLVENTATOMS,3,STOREDCONFS)
      double precision StoreSvnCentre(MAXSOLVENTS,3,STOREDCONFS)

      common /DSTORECOORDATA/ StoreZCoords,StoreZSolCoords,StoreSolCentre,
     . StoreSolGeomCentre,StoreSolDummies,StoreZSvnCoords,StoreSvnCentre
c   Keeping track of saved Born radii
      double precision StoreZSolContBorn(MAXSOLUTES,MAXSOLUTEATOMS,STOREDNRGS)
      double precision StoreZSolContBornF(MAXSOLUTES,MAXSOLUTEATOMS,STOREDNRGS)
      double precision StoreZSolContBornB(MAXSOLUTES,MAXSOLUTEATOMS,STOREDNRGS)
      double precision StoreZSolEffBorn(MAXSOLUTES,MAXSOLUTEATOMS,STOREDNRGS)
      double precision StoreZSolEffBornF(MAXSOLUTES,MAXSOLUTEATOMS,STOREDNRGS)
      double precision StoreZSolEffBornB(MAXSOLUTES,MAXSOLUTEATOMS,STOREDNRGS)
      double precision StoreZProtContBorn(MAXPROTEINS,MAXATOMS,STOREDNRGS)
      double precision StoreZProtContBornF(MAXPROTEINS,MAXATOMS,STOREDNRGS)
      double precision StoreZProtContBornB(MAXPROTEINS,MAXATOMS,STOREDNRGS)
      double precision StoreZProtEffBorn(MAXPROTEINS,MAXATOMS,STOREDNRGS)
      double precision StoreZProtEffBornF(MAXPROTEINS,MAXATOMS,STOREDNRGS)
      double precision StoreZProtEffBornB(MAXPROTEINS,MAXATOMS,STOREDNRGS)

      double precision StoreZSvnContBorn(MAXSOLVENTS,MAXSOLVENTATOMS,STOREDNRGS)
      double precision StoreZSvnContBornF(MAXSOLVENTS,MAXSOLVENTATOMS,STOREDNRGS)
      double precision StoreZSvnContBornB(MAXSOLVENTS,MAXSOLVENTATOMS,STOREDNRGS)
      double precision StoreZSvnEffBorn(MAXSOLVENTS,MAXSOLVENTATOMS,STOREDNRGS)
      double precision StoreZSvnEffBornF(MAXSOLVENTS,MAXSOLVENTATOMS,STOREDNRGS)
      double precision StoreZSvnEffBornB(MAXSOLVENTS,MAXSOLVENTATOMS,STOREDNRGS)
      double precision StoreOldZSolContBorn(MAXSOLUTES,MAXSOLUTEATOMS,STOREDNRGS)
      double precision StoreOldZSolContBornF(MAXSOLUTES,MAXSOLUTEATOMS,STOREDNRGS)
      double precision StoreOldZSolContBornB(MAXSOLUTES,MAXSOLUTEATOMS,STOREDNRGS)
      double precision StoreOldZSolEffBorn(MAXSOLUTES,MAXSOLUTEATOMS,STOREDNRGS)
      double precision StoreOldZSolEffBornF(MAXSOLUTES,MAXSOLUTEATOMS,STOREDNRGS)
      double precision StoreOldZSolEffBornB(MAXSOLUTES,MAXSOLUTEATOMS,STOREDNRGS)
      double precision StoreOldZProtContBorn(MAXPROTEINS,MAXATOMS,STOREDNRGS)
      double precision StoreOldZProtContBornF(MAXPROTEINS,MAXATOMS,STOREDNRGS)
      double precision StoreOldZProtContBornB(MAXPROTEINS,MAXATOMS,STOREDNRGS)
      double precision StoreOldZProtEffBorn(MAXPROTEINS,MAXATOMS,STOREDNRGS)
      double precision StoreOldZProtEffBornF(MAXPROTEINS,MAXATOMS,STOREDNRGS)
      double precision StoreOldZProtEffBornB(MAXPROTEINS,MAXATOMS,STOREDNRGS)
      double precision StoreOldZSvnContBorn(MAXSOLVENTS,MAXSOLVENTATOMS,STOREDNRGS)
      double precision StoreOldZSvnContBornF(MAXSOLVENTS,MAXSOLVENTATOMS,STOREDNRGS)
      double precision StoreOldZSvnContBornB(MAXSOLVENTS,MAXSOLVENTATOMS,STOREDNRGS)
      double precision StoreOldZSvnEffBorn(MAXSOLVENTS,MAXSOLVENTATOMS,STOREDNRGS)
      double precision StoreOldZSvnEffBornF(MAXSOLVENTS,MAXSOLVENTATOMS,STOREDNRGS)
      double precision StoreOldZSvnEffBornB(MAXSOLVENTS,MAXSOLVENTATOMS,STOREDNRGS)

      common /DSTOREBORNDATA/ StoreZSolContBorn,StoreZSolContBornF,
     . StoreZSolContBornB,StoreZProtContBorn,StoreZProtContBornF,
     . StoreZProtContBornB,StoreZSvnContBorn,StoreZSvnContBornF,
     . StoreZSvnContBornB,StoreOldZSolContBorn,StoreOldZSolContBornF,
     . StoreOldZSolContBornB,StoreOldZProtContBorn,StoreOldZProtContBornF,
     . StoreOldZProtContBornB,StoreOldZSvnContBorn,StoreOldZSvnContBornF,
     . StoreOldZSvnContBornB,StoreOldZSolEffBorn,StoreOldZSolEffBornF,
     . StoreOldZSolEffBornB,StoreOldZProtEffBorn,StoreOldZProtEffBornF,
     . StoreOldZProtEffBornB,StoreOldZSvnEffBorn,StoreOldZSvnEffBornF,
     . StoreOldZSvnEffBornB,StoreZSolEffBorn,StoreZSolEffBornF,
     . StoreZSolEffBornB,StoreZProtEffBorn,StoreZProtEffBornF,StoreZProtEffBornB,
     . StoreZSvnEffBorn,StoreZSvnEffBornF,StoreZSvnEffBornB
c
c     Variables keeping track of what has moved during a split move 
c     sequence
c
      integer MAXSPLITSOLMOVED
      parameter (MAXSPLITSOLMOVED=10)
      integer MAXSPLITRESMOVED
      parameter (MAXSPLITRESMOVED=50)
      integer MAXSPLITSVNMOVED
      parameter (MAXSPLITSVNMOVED=50)
      integer SplitSolMovedList(MAXSPLITSOLMOVED)
      integer nSolSplitMoved
      integer SplitResMovedList(MAXSPLITRESMOVED,2)
      integer nResSplitMoved
      integer SplitSvnMovedList(MAXSPLITSVNMOVED)
      integer nSvnSplitMoved

      common /ISPLITMOVE/ SplitSolMovedList,nSolSplitMoved,
     .  SplitResMovedList,nResSplitMoved,SplitSvnMovedList,
     .  nSvnSplitMoved
     
c     Whether or not the approximate potential has been activated
      logical QuickMove
c     Are we using a Distance Dependent Dielectric ?
      logical DDDon
c     Are we using constant born radii ? 
      logical ConstBorn

      common /LSPLITMOVE/ QuickMove,DDDon,ConstBorn



     
