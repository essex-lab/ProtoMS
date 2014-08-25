
c   This file holds all of the variables that are used to store the
c   energies

c       protein-protein energies
      double precision NewProProC(MAXPROTEINS,MAXPROTEINS)
      double precision NewProProLJ(MAXPROTEINS,MAXPROTEINS)
      double precision NewProProDist(MAXPROTEINS,MAXPROTEINS)
      double precision OldProProC(MAXPROTEINS,MAXPROTEINS)
      double precision OldProProLJ(MAXPROTEINS,MAXPROTEINS)
      double precision OldProProDist(MAXPROTEINS,MAXPROTEINS)
      double precision NewProProCF(MAXPROTEINS,MAXPROTEINS)
      double precision NewProProLJF(MAXPROTEINS,MAXPROTEINS)
      double precision OldProProCF(MAXPROTEINS,MAXPROTEINS)
      double precision OldProProLJF(MAXPROTEINS,MAXPROTEINS)
      double precision NewProProCB(MAXPROTEINS,MAXPROTEINS)
      double precision NewProProLJB(MAXPROTEINS,MAXPROTEINS)
      double precision OldProProCB(MAXPROTEINS,MAXPROTEINS)
      double precision OldProProLJB(MAXPROTEINS,MAXPROTEINS)
c       intramolecular protein energies
      double precision NewProBnd(MAXPROTEINS)
      double precision NewProAng(MAXPROTEINS)
      double precision NewProDih(MAXPROTEINS)
      double precision NewProUBy(MAXPROTEINS)
      double precision NewProNBC(MAXPROTEINS)
      double precision NewProNBLJ(MAXPROTEINS)
      double precision OldProBnd(MAXPROTEINS)
      double precision OldProAng(MAXPROTEINS)
      double precision OldProDih(MAXPROTEINS)
      double precision OldProUBy(MAXPROTEINS)
      double precision OldProNBC(MAXPROTEINS)
      double precision OldProNBLJ(MAXPROTEINS)
      double precision NewProBndF(MAXPROTEINS)
      double precision NewProAngF(MAXPROTEINS)
      double precision NewProDihF(MAXPROTEINS)
      double precision NewProUByF(MAXPROTEINS)
      double precision NewProNBCF(MAXPROTEINS)
      double precision NewProNBLJF(MAXPROTEINS)
      double precision OldProBndF(MAXPROTEINS)
      double precision OldProAngF(MAXPROTEINS)
      double precision OldProDihF(MAXPROTEINS)
      double precision OldProUByF(MAXPROTEINS)
      double precision OldProNBCF(MAXPROTEINS)
      double precision OldProNBLJF(MAXPROTEINS)
      double precision NewProBndB(MAXPROTEINS)
      double precision NewProAngB(MAXPROTEINS)
      double precision NewProDihB(MAXPROTEINS)
      double precision NewProUByB(MAXPROTEINS)
      double precision NewProNBCB(MAXPROTEINS)
      double precision NewProNBLJB(MAXPROTEINS)
      double precision OldProBndB(MAXPROTEINS)
      double precision OldProAngB(MAXPROTEINS)
      double precision OldProDihB(MAXPROTEINS)
      double precision OldProUByB(MAXPROTEINS)
      double precision OldProNBCB(MAXPROTEINS)
      double precision OldProNBLJB(MAXPROTEINS)
c       solute-solute energies
      double precision NewSolSolC(MAXSOLUTES,MAXSOLUTES)
      double precision NewSolSolLJ(MAXSOLUTES,MAXSOLUTES)
      double precision NewSolSolDist(MAXSOLUTES,MAXSOLUTES)
      double precision OldSolSolC(MAXSOLUTES,MAXSOLUTES)
      double precision OldSolSolLJ(MAXSOLUTES,MAXSOLUTES)
      double precision OldSolSolDist(MAXSOLUTES,MAXSOLUTES)
      double precision NewSolSolCF(MAXSOLUTES,MAXSOLUTES)
      double precision NewSolSolLJF(MAXSOLUTES,MAXSOLUTES)
      double precision OldSolSolCF(MAXSOLUTES,MAXSOLUTES)
      double precision OldSolSolLJF(MAXSOLUTES,MAXSOLUTES)
      double precision NewSolSolCB(MAXSOLUTES,MAXSOLUTES)
      double precision NewSolSolLJB(MAXSOLUTES,MAXSOLUTES)
      double precision OldSolSolCB(MAXSOLUTES,MAXSOLUTES)
      double precision OldSolSolLJB(MAXSOLUTES,MAXSOLUTES)
c       solute intramolecular energies
      double precision NewSolBnd(MAXSOLUTES)
      double precision NewSolAng(MAXSOLUTES)
      double precision NewSolDih(MAXSOLUTES)
      double precision NewSolUBy(MAXSOLUTES)
      double precision NewSolNBC(MAXSOLUTES)
      double precision NewSolNBLJ(MAXSOLUTES)
      double precision OldSolBnd(MAXSOLUTES)
      double precision OldSolAng(MAXSOLUTES)
      double precision OldSolDih(MAXSOLUTES)
      double precision OldSolUBy(MAXSOLUTES)
      double precision OldSolNBC(MAXSOLUTES)
      double precision OldSolNBLJ(MAXSOLUTES)
      double precision NewSolBndF(MAXSOLUTES)
      double precision NewSolAngF(MAXSOLUTES)
      double precision NewSolDihF(MAXSOLUTES)
      double precision NewSolUByF(MAXSOLUTES)
      double precision NewSolNBCF(MAXSOLUTES)
      double precision NewSolNBLJF(MAXSOLUTES)
      double precision OldSolBndF(MAXSOLUTES)
      double precision OldSolAngF(MAXSOLUTES)
      double precision OldSolDihF(MAXSOLUTES)
      double precision OldSolUByF(MAXSOLUTES)
      double precision OldSolNBCF(MAXSOLUTES)
      double precision OldSolNBLJF(MAXSOLUTES)
      double precision NewSolBndB(MAXSOLUTES)
      double precision NewSolAngB(MAXSOLUTES)
      double precision NewSolDihB(MAXSOLUTES)
      double precision NewSolUByB(MAXSOLUTES)
      double precision NewSolNBCB(MAXSOLUTES)
      double precision NewSolNBLJB(MAXSOLUTES)
      double precision OldSolBndB(MAXSOLUTES)
      double precision OldSolAngB(MAXSOLUTES)
      double precision OldSolDihB(MAXSOLUTES)
      double precision OldSolUByB(MAXSOLUTES)
      double precision OldSolNBCB(MAXSOLUTES)
      double precision OldSolNBLJB(MAXSOLUTES)
      double precision NewSoluteEnergy(MAXSOLUTES)
      double precision OldSoluteEnergy(MAXSOLUTES)       
      double precision NewSoluteEnergyF(MAXSOLUTES),NewSoluteEnergyB(MAXSOLUTES)
      double precision OldSoluteEnergyF(MAXSOLUTES),OldSoluteEnergyB(MAXSOLUTES)
c       solute-protein energies
      double precision NewProSolC(MAXPROTEINS,MAXSOLUTES)
      double precision NewProSolLJ(MAXPROTEINS,MAXSOLUTES)
      double precision NewProSolDist(MAXPROTEINS,MAXSOLUTES)
      double precision OldProSolC(MAXPROTEINS,MAXSOLUTES)
      double precision OldProSolLJ(MAXPROTEINS,MAXSOLUTES)
      double precision OldProSolDist(MAXPROTEINS,MAXSOLUTES)
      double precision NewProSolCF(MAXPROTEINS,MAXSOLUTES)
      double precision NewProSolLJF(MAXPROTEINS,MAXSOLUTES)
      double precision OldProSolCF(MAXPROTEINS,MAXSOLUTES)
      double precision OldProSolLJF(MAXPROTEINS,MAXSOLUTES)
      double precision NewProSolCB(MAXPROTEINS,MAXSOLUTES)
      double precision NewProSolLJB(MAXPROTEINS,MAXSOLUTES)
      double precision OldProSolCB(MAXPROTEINS,MAXSOLUTES)
      double precision OldProSolLJB(MAXPROTEINS,MAXSOLUTES)
c       solvent solvent energy
      double precision NewSvnSvnC
      double precision NewSvnSvnLJ
      double precision NewSvnSvnDist
      double precision OldSvnSvnC
      double precision OldSvnSvnLJ
      double precision OldSvnSvnDist
      double precision NewSvnSvnCF
      double precision NewSvnSvnLJF
      double precision OldSvnSvnCF
      double precision OldSvnSvnLJF
      double precision NewSvnSvnCB
      double precision NewSvnSvnLJB
      double precision OldSvnSvnCB
      double precision OldSvnSvnLJB



c       solute-solvent energy
      double precision NewSolSvnC(MAXSOLUTES)
      double precision NewSolSvnLJ(MAXSOLUTES)
      double precision NewSolSvnDist(MAXSOLUTES)
      double precision OldSolSvnC(MAXSOLUTES)
      double precision OldSolSvnLJ(MAXSOLUTES)
      double precision OldSolSvnDist(MAXSOLUTES)
      double precision NewSolSvnCF(MAXSOLUTES)
      double precision NewSolSvnLJF(MAXSOLUTES)
      double precision OldSolSvnCF(MAXSOLUTES)
      double precision OldSolSvnLJF(MAXSOLUTES)
      double precision NewSolSvnCB(MAXSOLUTES)
      double precision NewSolSvnLJB(MAXSOLUTES)
      double precision OldSolSvnCB(MAXSOLUTES)
      double precision OldSolSvnLJB(MAXSOLUTES)
c      solute-gcs energy
      double precision NewSolGCSC(MAXSOLUTES)
      double precision NewSolGCSLJ(MAXSOLUTES)
      double precision OldSolGCSC(MAXSOLUTES)
      double precision OldSolGCSLJ(MAXSOLUTES)
      double precision NewSolGCSCF(MAXSOLUTES)
      double precision NewSolGCSLJF(MAXSOLUTES)
      double precision OldSolGCSCF(MAXSOLUTES)
      double precision OldSolGCSLJF(MAXSOLUTES)
      double precision NewSolGCSCB(MAXSOLUTES)
      double precision NewSolGCSLJB(MAXSOLUTES)
      double precision OldSolGCSCB(MAXSOLUTES)
      double precision OldSolGCSLJB(MAXSOLUTES)
c       gcsolute-solvent energies
      double precision NewGCSSvnC(MAXSOLVENTS)
      double precision NewGCSSvnLJ(MAXSOLVENTS)
      double precision OldGCSSvnC(MAXSOLVENTS)
      double precision OldGCSSvnLJ(MAXSOLVENTS)
      double precision NewGCSSvnCF(MAXSOLVENTS)
      double precision NewGCSSvnLJF(MAXSOLVENTS)
      double precision OldGCSSvnCF(MAXSOLVENTS)
      double precision OldGCSSvnLJF(MAXSOLVENTS)
      double precision NewGCSSvnCB(MAXSOLVENTS)
      double precision NewGCSSvnLJB(MAXSOLVENTS)
      double precision OldGCSSvnCB(MAXSOLVENTS)
      double precision OldGCSSvnLJB(MAXSOLVENTS)
c       gcs-solute energy (ind)
      double precision NewIndGCSSolC(MAXSOLVENTS)
      double precision OldIndGCSSolC(MAXSOLVENTS)
      double precision NewIndGCSSolCF(MAXSOLVENTS)
      double precision OldIndGCSSolCF(MAXSOLVENTS)
      double precision NewIndGCSSolCB(MAXSOLVENTS)
      double precision OldIndGCSSolCB(MAXSOLVENTS)
      double precision NewIndGCSSolLJ(MAXSOLVENTS)
      double precision OldIndGCSSolLJ(MAXSOLVENTS)
      double precision NewIndGCSSolLJB(MAXSOLVENTS)
      double precision OldIndGCSSolLJB(MAXSOLVENTS)
      double precision NewIndGCSSolLJF(MAXSOLVENTS)
      double precision OldIndGCSSolLJF(MAXSOLVENTS)
c       protein-solvent energy
      double precision NewProSvnC(MAXPROTEINS)
      double precision NewProSvnLJ(MAXPROTEINS)
      double precision NewProSvnDist(MAXPROTEINS)
      double precision OldProSvnC(MAXPROTEINS)
      double precision OldProSvnLJ(MAXPROTEINS)
      double precision OldProSvnDist(MAXPROTEINS)
      double precision NewProSvnCF(MAXPROTEINS)
      double precision NewProSvnLJF(MAXPROTEINS)
      double precision OldProSvnCF(MAXPROTEINS)
      double precision OldProSvnLJF(MAXPROTEINS)
      double precision NewProSvnCB(MAXPROTEINS)
      double precision NewProSvnLJB(MAXPROTEINS)
      double precision OldProSvnCB(MAXPROTEINS)
      double precision OldProSvnLJB(MAXPROTEINS)
c      protein-gcs energy
      double precision NewProGCSC(MAXPROTEINS)
      double precision NewProGCSLJ(MAXPROTEINS)
      double precision NewProGCSDist(MAXPROTEINS)
      double precision OldProGCSC(MAXPROTEINS)
      double precision OldProGCSLJ(MAXPROTEINS)
      double precision OldProGCSDist(MAXPROTEINS)
      double precision NewProGCSCF(MAXPROTEINS)
      double precision NewProGCSLJF(MAXPROTEINS)
      double precision OldProGCSCF(MAXPROTEINS)
      double precision OldProGCSLJF(MAXPROTEINS)
      double precision NewProGCSCB(MAXPROTEINS)
      double precision NewProGCSLJB(MAXPROTEINS)
      double precision OldProGCSCB(MAXPROTEINS)
      double precision OldProGCSLJB(MAXPROTEINS)
c       protein residue energies
      double precision oldResProC(MAXPROTEINS),oldResProL(MAXPROTEINS)
      double precision oldResSolC(MAXSOLUTES),oldResSolL(MAXSOLUTES)
      double precision oldResProCF(MAXPROTEINS),oldResProLF(MAXPROTEINS)
      double precision oldResSolCF(MAXSOLUTES),oldResSolLF(MAXSOLUTES)
      double precision oldResProCB(MAXPROTEINS),oldResProLB(MAXPROTEINS)
      double precision oldResSolCB(MAXSOLUTES),oldResSolLB(MAXSOLUTES)
      double precision oldResProDist(MAXPROTEINS),newResProDist(MAXPROTEINS)
      double precision oldResSolDist(MAXSOLUTES),newResSolDist(MAXSOLUTES)
      double precision oldResSvnDist,newResSvnDist
      double precision oldResGCSDist,newResGCSDist
      double precision oldResSvnC,oldResSvnL
      double precision oldResBnd,oldResAng,oldResDih,oldResUBy
      double precision oldResNBC,oldResNBLJ
      double precision newResProC(MAXPROTEINS),newResProL(MAXPROTEINS)
      double precision newResSolC(MAXSOLUTES),newResSolL(MAXSOLUTES)
      double precision newResSvnC,newResSvnL
      double precision newResBnd,newResAng,newResDih,newResUBy
      double precision newResNBC,newResNBLJ
      double precision newResNrg,oldResNrg
      double precision newResSvnCF,newResSvnLF
      double precision newResBndF,newResAngF,newResDihF,newResUByF
      double precision newResNBCF,newResNBLJF
      double precision oldResSvnCF,oldResSvnLF
      double precision oldResBndF,oldResAngF,oldResDihF,oldResUByF
      double precision oldResNBCF,oldResNBLJF
      double precision newResProCF(MAXPROTEINS),newResProLF(MAXPROTEINS)
      double precision newResSolCF(MAXSOLUTES),newResSolLF(MAXSOLUTES)
      double precision newResSvnCB,newResSvnLB
      double precision newResBndB,newResAngB,newResDihB,newResUByB
      double precision newResNBCB,newResNBLJB
      double precision oldResSvnCB,oldResSvnLB
      double precision oldResBndB,oldResAngB,oldResDihB,oldResUByB
      double precision oldResNBCB,oldResNBLJB
      double precision newResProCB(MAXPROTEINS),newResProLB(MAXPROTEINS)
      double precision newResSolCB(MAXSOLUTES),newResSolLB(MAXSOLUTES)
      double precision newResGCSC,newResGCSL
      double precision newResGCSCF,newResGCSLF
      double precision newResGCSCB,newResGCSLB
      double precision oldResGCSC,oldResGCSL
      double precision oldResGCSCF,oldResGCSLF
      double precision oldResGCSCB,oldResGCSLB
c       individual solvent - other energies
      double precision oldIndSvnGCSC(MAXSOLVENTS),newIndSvnGCSC(MAXSOLVENTS)
      double precision oldIndSvnGCSCF(MAXSOLVENTS),newIndSvnGCSCF(MAXSOLVENTS)
      double precision oldIndSvnGCSCB(MAXSOLVENTS),newIndSvnGCSCB(MAXSOLVENTS)
      double precision oldIndSvnGCSLJ(MAXSOLVENTS),newIndSvnGCSLJ(MAXSOLVENTS)
      double precision oldIndSvnGCSLJF(MAXSOLVENTS),newIndSvnGCSLJF(MAXSOLVENTS)
      double precision oldIndSvnGCSLJB(MAXSOLVENTS),newIndSvnGCSLJB(MAXSOLVENTS)
      double precision oldIndSvnSvnC,newIndSvnSvnC
      double precision oldIndSvnProC(MAXPROTEINS),newIndSvnProC(MAXPROTEINS)
      double precision oldIndSvnSolC(MAXSOLUTES),newIndSvnSolC(MAXSOLUTES)
      double precision oldIndSvnC,newIndSvnC
      double precision oldIndSvnSvnLJ,newIndSvnSvnLJ
      double precision oldIndSvnProLJ(MAXPROTEINS),newIndSvnProLJ(MAXPROTEINS)
      double precision oldIndSvnSolLJ(MAXSOLUTES),newIndSvnSolLJ(MAXSOLUTES)
      double precision oldIndSvnLJ,newIndSvnLJ
      double precision oldIndSvnSvnCF,newIndSvnSvnCF
      double precision oldIndSvnProCF(MAXPROTEINS),newIndSvnProCF(MAXPROTEINS)
      double precision oldIndSvnSolCF(MAXSOLUTES),newIndSvnSolCF(MAXSOLUTES)
      double precision oldIndSvnCF,newIndSvnCF
      double precision oldIndSvnSvnLJF,newIndSvnSvnLJF
      double precision oldIndSvnProLJF(MAXPROTEINS),newIndSvnProLJF(MAXPROTEINS)
      double precision oldIndSvnSolLJF(MAXSOLUTES),newIndSvnSolLJF(MAXSOLUTES)
      double precision oldIndSvnLJF,newIndSvnLJF
      double precision oldIndSvnSvnCB,newIndSvnSvnCB
      double precision oldIndSvnProCB(MAXPROTEINS),newIndSvnProCB(MAXPROTEINS)
      double precision oldIndSvnSolCB(MAXSOLUTES),newIndSvnSolCB(MAXSOLUTES)
      double precision oldIndSvnCB,newIndSvnCB
      double precision oldIndSvnSvnLJB,newIndSvnSvnLJB
      double precision oldIndSvnProLJB(MAXPROTEINS),newIndSvnProLJB(MAXPROTEINS)
      double precision oldIndSvnSolLJB(MAXSOLUTES),newIndSvnSolLJB(MAXSOLUTES)
      double precision oldIndSvnLJB,newIndSvnLJB
c     old and new residue energies
      double precision oldResNrgF,oldResNrgB
      double precision newResNrgF,newResNrgB
c       total energies!
      double precision NewTotalEnergy
      double precision OldTotalEnergy      
      double precision OldTotalEnergyF,OldTotalEnergyB
      double precision NewTotalEnergyF,NewTotalEnergyB
c      GC-solute / GC-solute
      double precision newGCSGCSC,oldGCSGCSC
      double precision newGCSGCSCF,oldGCSGCSCF
      double precision newGCSGCSCB,oldGCSGCSCB
      double precision newGCSGCSLJ,oldGCSGCSLJ
      double precision newGCSGCSLJB,oldGCSGCSLJB
      double precision newGCSGCSLJF,oldGCSGCSLJF
      
      double precision newIndGCSGCSC,oldIndGCSGCSC
      double precision newIndGCSGCSCF,oldIndGCSGCSCF
      double precision newIndGCSGCSCB,oldIndGCSGCSCB
      double precision newIndGCSGCSLJ,oldIndGCSGCSLJ
      double precision newIndGCSGCSLJB,oldIndGCSGCSLJB
      double precision newIndGCSGCSLJF,oldIndGCSGCSLJF
       
c        GC-solute / Protein
      double precision newIndGCSProC(MAXPROTEINS),oldIndGCSProC(MAXPROTEINS)
      double precision newIndGCSProCF(MAXPROTEINS),oldIndGCSProCF(MAXPROTEINS)
      double precision newIndGCSProCB(MAXPROTEINS),oldIndGCSProCB(MAXPROTEINS)
      double precision newIndGCSProLJ(MAXPROTEINS),oldIndGCSProLJ(MAXPROTEINS)
      double precision newIndGCSProLJB(MAXPROTEINS),oldIndGCSProLJB(MAXPROTEINS)
      double precision newIndGCSProLJF(MAXPROTEINS),oldIndGCSProLJF(MAXPROTEINS)
c      GC-solute / Solvent
      double precision newIndGCSSvnC(MAXSOLVENTS),oldIndGCSSvnC(MAXSOLVENTS)
      double precision newIndGCSSvnCF(MAXSOLVENTS),oldIndGCSSvnCF(MAXSOLVENTS)
      double precision newIndGCSSvnCB(MAXSOLVENTS),oldIndGCSSvnCB(MAXSOLVENTS)
      double precision newIndGCSSvnLJ(MAXSOLVENTS),oldIndGCSSvnLJ(MAXSOLVENTS)
      double precision newIndGCSSvnLJB(MAXSOLVENTS),oldIndGCSSvnLJB(MAXSOLVENTS)
      double precision newIndGCSSvnLJF(MAXSOLVENTS),oldIndGCSSvnLJF(MAXSOLVENTS)
      double precision oldGCSSvnDist(MAXSOLVENTS),newGCSSvnDist(MAXSOLVENTS)
c      Overall
      double precision oldIndGCSC,newIndGCSC
      double precision oldIndGCSLJ,newIndGCSLJ
      double precision oldIndGCSCF,newIndGCSCF
      double precision oldIndGCSLJF,newIndGCSLJF
      double precision oldIndGCSCB,newIndGCSCB
      double precision oldIndGCSLJB,newIndGCSLJB
      
c     Energies for MBAR, BAR or TI
      double precision MBAREnergy(MAXLAMBDAS),AnalyticalGradient,NumericalGradient
  

c     the common blocks...
      common /DENERGY1/ NewProProC,NewProProLJ,NewProBnd,NewProAng,
     .  NewProDih,NewProUBy,NewProNBC,NewProNBLJ,NewSolSolC,NewSolSolLJ,
     .  NewSolBnd,NewSolAng,NewSolDih,NewSolUBy,NewSolNBC,NewSolNBLJ,
     .  NewProSolC,NewProSolLJ,NewSvnSvnC,NewSvnSvnLJ,
     .  NewSolSvnC,NewSolSvnLJ,NewProSvnC,NewProSvnLJ,NewTotalEnergy,
     .  NewProProDist,NewProSolDist,NewProSvnDist,
     .  NewSolSolDist,NewSolSvnDist,NewSvnSvnDist,
     .  OldProProC,OldProProLJ,OldProBnd,OldProAng,
     .  OldProDih,OldProUBy,OldProNBC,OldProNBLJ,OldSolSolC,OldSolSolLJ,
     .  OldSolBnd,OldSolAng,OldSolDih,OldSolUBy,OldSolNBC,OldSolNBLJ,
     .  OldProSolC,OldProSolLJ,OldSvnSvnC,OldSvnSvnLJ,
     .  OldSolSvnC,OldSolSvnLJ,OldProSvnC,OldProSvnLJ,OldTotalEnergy,
     .  OldProProDist,OldProSolDist,OldProSvnDist,NewSoluteEnergy,
     .  OldSolSolDist,OldSolSvnDist,OldSvnSvnDist,OldSoluteEnergy,
     .  OldTotalEnergyF,OldTotalEnergyB,NewTotalEnergyF,NewTotalEnergyB,
     .  NewProGCSC,NewProGCSLJ,NewProGCSCF,NewProGCSLJF,NewProGCSCB,NewProGCSLJB,
     .  OldProGCSC,OldProGCSLJ,OldProGCSCF,OldProGCSCB,OldProGCSLJB,OldProGCSLJF,
     .  NewProGCSDist,OldProGCSDist,newIndSvnGCSLJF,newIndSvnGCSLJB,oldIndSvnGCSLJF,oldIndSvnGCSLJB,
     .  newIndSvnGCSLJ,oldIndSvnGCSLJ,oldIndSvnGCSC,newIndSvnGCSC,oldIndSvnGCSCF,newIndSvnGCSCF,
     .  oldIndSvnGCSCB,newIndSvnGCSCB,NewGCSSvnC,OldGCSSvnC,NewGCSSvnCF,OldGCSSvnCF,NewGCSSvnCB,OldGCSSvnCB,
     .  NewGCSSvnLJ,oldGCSSvnLJ,NewGCSSvnLJF,oldGCSSvnLJF,NewGCSSvnLJB,oldGCSSvnLJB,
     .  OldSolGCSC,NewSolGCSC,OldSolGCSCF,NewSolGCSCF,OldSolGCSCB,NewSolGCSCB,
     .  OldSolGCSLJ,NewSolGCSLJ,OldSolGCSLJF,NewSolGCSLJF,OldSolGCSLJB,NewSolGCSLJB

     
 
     
      common /DENERGY2/ oldResProC,oldResProL,oldResSolC,oldResSolL,
     .                      oldResSvnC,oldResSvnL,oldResBnd,oldResAng,oldResDih,
     . oldResUBy,oldResNBC,oldResNBLJ,newResProC,newResProL,
     . newResSolC,newResSolL,newResSvnC,newResSvnL,newResBnd,newResAng,newResDih,
     . newResUBy,newResNBC,newResNBLJ,oldResNrg,newResNrg,oldresprodist,newresprodist,
     . oldressoldist,newressoldist,oldressvndist,newressvndist,
     . oldIndSvnSvnC,newIndSvnSvnC,oldIndSvnSvnLJ,newIndSvnSvnLJ,
     . oldIndSvnProC,newIndSvnProC,oldIndSvnProLJ,newIndSvnProLJ,
     . oldIndSvnSolC,newIndSvnSolC,oldIndSvnSolLJ,newIndSvnSolLJ,
     . oldIndSvnC,oldIndSvnLJ,newIndSvnC,newIndSvnLJ,
     . oldResNrgF,oldResNrgB,newResNrgF,newResNrgB,
     . oldIndSvnLJF,oldIndSvnLJB,oldIndSvnCF,oldIndSvnCB,
     . newIndSvnLJF,newIndSvnLJB,newIndSvnCF,newIndSvnCB,oldResGCSDist,newResGCSDist,
     . newSoluteEnergyF,newSoluteEnergyB,oldSoluteEnergyF,oldSoluteEnergyB,
     . newResGCSC,newResGCSL,newResGCSCF,newResGCSLF,newResGCSCB,newResGCSLB,
     . oldResGCSC,oldResGCSL,oldResGCSCF,oldResGCSLF,oldResGCSCB,oldResGCSLB,
     . NewIndGCSSolC,OldIndGCSSolC,NewIndGCSSolCF,OldIndGCSSolCF,NewIndGCSSolCB,OldIndGCSSolCB,
     . NewIndGCSSolLJ,OldIndGCSSolLJ,NewIndGCSSolLJF,OldIndGCSSolLJF,NewIndGCSSolLJB,OldIndGCSSolLJB

      common /DENERGY3/ newindsvnsolljb,oldindsvnsolljb,newindsvnproljb,
     . oldindsvnproljb,newindsvnsvnljb,oldindsvnsvnljb,newindsvnsolcb,
     . oldindsvnsolcb,newindsvnprocb,oldindsvnprocb,newindsvnsvncb,
     . oldindsvnsvncb,newindsvnsolljf,oldindsvnsolljf,newindsvnproljf,
     . oldindsvnproljf,newindsvnsvnljf,oldindsvnsvnljf,newindsvnsolcf,
     . oldindsvnsolcf,newindsvnprocf,oldindsvnprocf,newindsvnsvncf
      

      common /DENERGY4/ oldindsvnsvncf,newressollb,newressolcb,newresprolb,newresprocb,
     . oldresnbljb,oldresnbcb,oldresdihb,oldresangb,oldresbndb,oldressvnlb,
     . oldresubyb,oldresubyf,newresubyb,newresubyf,
     . oldressvncb,newresnbljb,newresnbcb,newresdihb,newresangb,newresbndb,
     . newressvnlb,newressvncb,newressollf,newressolcf,newresprolf,newresprocf,
     . oldresnbljf,oldresnbcf,oldresdihf,oldresangf,oldresbndf,oldressvnlf,
     . oldressvncf,newresnbljf,newresnbcf,newresdihf,newresangf,newresbndf,
     . newressvnlf,newressvncf,oldressollb,oldressolcb,oldresprolb,oldresprocb,
     . oldressollf,oldressolcf,oldresprolf,oldresprocf,newprosvnljb,newprosvncb,
     . newprosvnljf,newprosvncf,newsolsvnljb,newsolsvncb,newsolsvnljf,newsolsvncf,
     . newsvnsvnljb,newsvnsvncb,newsvnsvnljf,newsvnsvncf,newprosolljb,newprosolcb,
     . newprosolljf,newprosolcf,newsolnbljb,newsolnbcb,newsoldihb,newsolangb,
     . newsolbndb,newsolnbljf,newsolnbcf,newsoldihf,newsolangf,newsolbndf,
     . newsolubyf,newsolubyb,oldsolubyf,oldsolubyb,
     . newsolsolljb,newsolsolcb,newsolsolljf,newsolsolcf,newpronbljb,newpronbcb,
     . newprodihb,newproangb,newprobndb,newpronbljf,newpronbcf,newprodihf,
     . newproangf,newprobndf,newproproljb,newproprocb,newproproljf,newproprocf

      common /DENERGY5/ oldprosvnljb,oldprosvncb,oldprosvnljf,oldprosvncf,oldsolsvnljb,
     .  oldsolsvncb,oldsolsvnljf,oldsolsvncf,oldsvnsvnljb,oldsvnsvncb,
     .  oldsvnsvnljf,oldsvnsvncf,oldprosolljb,oldprosolcb,oldprosolljf,
     .  oldprosolcf,oldsolnbljb,oldsolnbcb,oldsoldihb,oldsolangb,
     .  oldsolbndb,oldsolnbljf,oldsolnbcf,oldsoldihf,oldsolangf,
     .  oldsolbndf,oldsolsolljb,oldsolsolcb,oldsolsolljf,oldsolsolcf,
     .  oldpronbljb,oldpronbcb,oldprodihb,oldproangb,oldprobndb,
     .  oldpronbljf,oldpronbcf,oldprodihf,oldproangf,oldprobndf,
     .  oldproproljb,oldproprocb,oldproproljf,oldproprocf,
     .  oldproubyf,oldproubyb,newproubyb,newproubyf

      common /DENERGY6/ oldIndGCSC,newIndGCSC,oldIndGCSLJ,newIndGCSLJ,oldIndGCSCF,newIndGCSCF,
     .  oldIndGCSLJF,newIndGCSLJF,oldIndGCSCB,newIndGCSCB,oldIndGCSLJB,newIndGCSLJB,
     .  newIndGCSSvnC,oldIndGCSSvnC,newIndGCSSvnCF,oldIndGCSSvnCF,newIndGCSSvnCB,oldIndGCSSvnCB,
     .  newIndGCSSvnLJ,oldIndGCSSvnLJ,newIndGCSSvnLJB,oldIndGCSSvnLJB,newIndGCSSvnLJF,oldIndGCSSvnLJF,
     .  newIndGCSProC,oldIndGCSProC,newIndGCSProCF,oldIndGCSProCF,newIndGCSProCB,oldIndGCSProCB,
     .  newIndGCSProLJ,oldIndGCSProLJ,newIndGCSProLJB,oldIndGCSProLJB,newIndGCSProLJF,oldIndGCSProLJF,
     .  newGCSGCSC,oldGCSGCSC,newGCSGCSCF,oldGCSGCSCF,newGCSGCSCB,oldGCSGCSCB,newIndGCSGCSLJB,oldIndGCSGCSLJB,
     .  newGCSGCSLJ,oldGCSGCSLJ,newGCSGCSLJB,oldGCSGCSLJB,newGCSGCSLJF,oldGCSGCSLJF,
     .  newGCSSvnDist,oldGCSSvnDist, newIndGCSGCSC,oldIndGCSGCSC, newIndGCSGCSLJ,oldIndGCSGCSLJ,
     .  newIndGCSGCSCF,oldIndGCSGCSCF, newIndGCSGCSCB,oldIndGCSGCSCB, newIndGCSGCSLJF,oldIndGCSGCSLJF
     
      common /DENERGY7/ MBAREnergy,AnalyticalGradient,NumericalGradient

      save /DENERGY1/
      save /DENERGY2/
      save /DENERGY3/
      save /DENERGY4/
      save /DENERGY5/
      save /DENERGY6/
      
