      subroutine revertEnergies(pot)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'solvents.inc'
      include 'flags.inc'
      include 'energies.inc'
      include 'printer.inc'
      include 'quickpot.inc'
      include 'constants.inc'
      include 'extraenergies.inc'
c########################################
c
c     Reverts the energies of a configuration
c     
c
c     (c) Julien Michel April 2005
c#########################################
      integer i,j,pot
c      write(printstring,*) 'REVERTING NERGIES FOR POT ',pot
c      call printLine(DEBUG,printstring)

      if (NProteins .gt. 0) then
        do i=1,NProteins
           OldProBnd(i) = StoredProBnd(i,pot)
           OldProAng(i) = StoredProAng(i,pot) 
           OldProDih(i) = StoredProDih(i,pot) 
           OldProUBy(i) = StoredProUBy(i,pot) 
           OldProNBC(i) = StoredProNBC(i,pot) 
           OldProNBLJ(i) = StoredProNBLJ(i,pot) 
           OldProBndF(i) = StoredProBndF(i,pot)
           OldProAngF(i) = StoredProAngF(i,pot)
           OldProDihF(i) = StoredProDihF(i,pot)
           OldProUByF(i) = StoredProUByF(i,pot)
           OldProNBCF(i) = StoredProNBCF(i,pot) 
           OldProNBLJF(i) = StoredProNBLJF(i,pot) 
           OldProBndB(i) = StoredProBndB(i,pot)
           OldProAngB(i) = StoredProAngB(i,pot) 
           OldProDihB(i) = StoredProDihB(i,pot) 
           OldProUByB(i) = StoredProUByB(i,pot) 
           OldProNBCB(i) = StoredProNBCB(i,pot) 
           OldProNBLJB(i) = StoredProNBLJB(i,pot) 
        enddo
        if (NProteins.gt.1) then
          do i=1,NProteins-1
            do j=i+1,NProteins
              OldProProC(i,j) = StoredProProC(i,j,pot) 
              OldProProLJ(i,j) = StoredProProLJ(i,j,pot)
              OldProProC(j,i) = StoredProProC(j,i,pot) 
              OldProProLJ(j,i) = StoredProProLJ(j,i,pot) 
              OldProProCF(i,j) = StoredProProCF(i,j,pot) 
              OldProProLJF(i,j) = StoredProProLJF(i,j,pot) 
              OldProProCF(j,i) = StoredProProCF(j,i,pot)
              OldProProLJF(j,i) = StoredProProLJF(j,i,pot)
              OldProProCB(i,j) = StoredProProCB(i,j,pot) 
              OldProProLJB(i,j) = StoredProProLJB(i,j,pot)
              OldProProCB(j,i) = StoredProProCB(j,i,pot)
              OldProProLJB(j,i) = StoredProProLJB(j,i,pot)
              OldProProDist(j,i) = StoredProProDist(j,i,pot) 
              OldProProDist(i,j) = StoredProProDist(i,j,pot)
            enddo
          enddo
        endif
      endif

      if (NSolutes.gt.0) then
        do i=1,NSolutes
          if (.not.IgnoreSolute(i)) then
            OldSolBnd(i) = StoredSolBnd(i,pot)
            OldSolAng(i) = StoredSolAng(i,pot)
            OldSolDih(i) = StoredSolDih(i,pot) 
            OldSolUBy(i) = StoredSolUBy(i,pot) 
            OldSolNBC(i) = StoredSolNBC(i,pot) 
            OldSolNBLJ(i) = StoredSolNBLJ(i,pot)
            OldSolBndF(i) = StoredSolBndF(i,pot) 
            OldSolAngF(i) = StoredSolAngF(i,pot)
            OldSolDihF(i) = StoredSolDihF(i,pot)
            OldSolUByF(i) = StoredSolUByF(i,pot) 
            OldSolNBCF(i) = StoredSolNBCF(i,pot) 
            OldSolNBLJF(i) = StoredSolNBLJ(i,pot)
            OldSolBndB(i) = StoredSolBndB(i,pot)
            OldSolAngB(i) = StoredSolAngB(i,pot)
            OldSolDihB(i) = StoredSolDihB(i,pot) 
            OldSolUByB(i) = StoredSolUByB(i,pot)
            OldSolNBCB(i) = StoredSolNBCB(i,pot) 
            OldSolNBLJB(i) = StoredSolNBLJB(i,pot) 
          endif
        enddo
        
        if (NSolutes.gt.1) then
          do i=1,NSolutes-1
            if (.not.IgnoreSolute(i)) then
              do j=i+1,NSolutes
                if (.not.IgnoreSolute(j)) then 
                  OldSolSolC(j,i) = StoredSolSolC(j,i,pot)
                  OldSolSolLJ(j,i) = StoredSolSolLJ(j,i,pot)
                  OldSolSolC(i,j) = StoredSolSolC(i,j,pot)
                  OldSolSolLJ(i,j) = StoredSolSolLJ(i,j,pot)
                  OldSolSolCF(j,i) = StoredSolSolCF(j,i,pot)
                  OldSolSolLJF(j,i) = StoredSolSolLJF(j,i,pot) 
                  OldSolSolCF(i,j) = StoredSolSolCF(i,j,pot)
                  OldSolSolLJF(i,j) = StoredSolSolLJF(i,j,pot)
                  OldSolSolCB(j,i) = StoredSolSolCB(j,i,pot) 
                  OldSolSolLJB(j,i) = StoredSolSolLJB(j,i,pot)
                  OldSolSolCB(i,j) = StoredSolSolCB(i,j,pot) 
                  OldSolSolLJB(i,j) = StoredSolSolLJB(i,j,pot) 
                  OldSolSolDist(j,i) = StoredSolSolDist(j,i,pot) 
                  OldSolSolDist(i,j) = StoredSolSolDist(i,j,pot) 
                endif
              enddo
            endif
          enddo
        endif

        if (NProteins.gt.0) then
          do j=1,NSolutes   
            if (.not.IgnoreSolute(j)) then
              do i=1,NProteins
                OldProSolC(i,j) = StoredProSolC(i,j,pot) 
                OldProSolLJ(i,j) = StoredProSolLJ(i,j,pot) 
                OldProSolCF(i,j) = StoredProSolCF(i,j,pot) 
                OldProSolLJF(i,j) = StoredProSolLJF(i,j,pot)
                OldProSolCB(i,j) = StoredProSolCB(i,j,pot) 
                OldProSolLJB(i,j) = StoredProSolLJB(i,j,pot) 
                OldProSolDist(i,j) = StoredProSolDist(i,j,pot) 
              enddo
            endif
          enddo
        endif

      endif

      if (NSolvents .gt. 0) then
        OldSvnSvnC = StoredSvnSvnC(pot) 
        OldSvnSvnLJ = StoredSvnSvnLJ(pot) 
        OldSvnSvnCF = StoredSvnSvnCF(pot) 
        OldSvnSvnLJF = StoredSvnSvnLJF(pot)
        OldSvnSvnCB = StoredSvnSvnCB(pot) 
        OldSvnSvnLJB = StoredSvnSvnLJB(pot)
        OldSvnSvnDist = StoredSvnSvnDist(pot)
        if (NProteins.gt.0) then
          do i=1,NProteins
            OldProSvnC(i) = StoredProSvnC(i,pot)
            OldProSvnLJ(i) = StoredProSvnLJ(i,pot) 
            OldProSvnCF(i) = StoredProSvnCF(i,pot) 
            OldProSvnLJF(i) = StoredProSvnLJF(i,pot) 
            OldProSvnCB(i)  = StoredProSvnCB(i,pot) 
            OldProSvnLJB(i) = StoredProSvnLJB(i,pot) 
            OldProSvnDist(i) = StoredProSvnDist(i,pot)
          enddo
        endif
        if (NSolutes.gt.0) then
          do i=1,NSolutes
            if (.not.IgnoreSolute(i)) then
              OldSolSvnC(i) = StoredSolSvnC(i,pot) 
              OldSolSvnLJ(i) = StoredSolSvnLJ(i,pot) 
              OldSolSvnCF(i) = StoredSolSvnCF(i,pot) 
              OldSolSvnLJF(i) = StoredSolSvnLJF(i,pot) 
              OldSolSvnCB(i) = StoredSolSvnCB(i,pot) 
              OldSolSvnLJB(i) = StoredSolSvnLJB(i,pot) 
              OldSolSvnDist(i) = StoredSolSvnDist(i,pot) 
              if (i.eq.SvnPrefSol .and. PrefSample) then
c               save the solvent probabilities...
                do j=1,NSolvents
                  SvnProb(j)  = StoredSvnProb(j,pot)
                enddo
                SvnSumProb = StoredSumProb(pot) 
                SvnMaxProb = StoredMaxProb(pot)
              endif
            endif
          enddo
        endif
      endif

      if (SolventCap) then 
        oldSvnCapEnergyPart = StoredSolventCapPart(pot) 
        oldSvnCapEnergy = StoredSolventCap(pot) 
      endif

      if (NHarmonicRestraints .gt. 0) then 
        oldHarmEnergyPart = StoredHarmEnergyPart(pot)
        oldHarmEnergyPartB = StoredHarmEnergyPartB(pot) 
        oldHarmEnergyPartF = StoredHarmEnergyPartF(pot) 
        oldHarmEnergy = StoredHarmEnergy(pot)
        oldHarmEnergyB = StoredHarmEnergyB(pot) 
        oldHarmEnergyF = StoredHarmEnergyF(pot) 
      endif

      if (NFlatBottomRestraints .gt. 0) then
        oldFBotEnergyPart = StoredFBotEnergyPart(pot) 
        oldFBotEnergyPartB = StoredFBotEnergyPartB(pot) 
        oldFBotEnergyPartF = StoredFBotEnergyPartF(pot) 
        oldFBotEnergy = StoredFBotEnergy(pot) 
        oldFBotEnergyB = StoredFBotEnergyB(pot) 
        oldFBotEnergyF = StoredFBotEnergyF(pot) 
      endif

      if (doSA) then
        oldTotalSurface = StoredTotalSurface(pot) 
        oldTotalSurfaceF = StoredTotalSurfaceF(pot) 
        oldTotalSurfaceB = StoredTotalSurfaceB(pot)
        oldSAEnergy = StoredSAEnergy(pot) 
        oldSAEnergyF = StoredSAEnergyF(pot) 
        oldSAEnergyB = StoredSAEnergyB(pot) 
      endif

      if (doGB) then
c        write(printstring,*) 'Revert GB ',StoredGBEnergy(pot),' pot ',pot
c        call printLine(DEBUG,printstring)
        OldGBEnergy = StoredGBEnergy(pot) 
        OldGBEnergyF = StoredGBEnergyF(pot) 
        OldGBEnergyB = StoredGBEnergyB(pot) 
      endif

      oldextraEnergyPart = StoredExtraEnergyPart(pot) 
      oldextraEnergyPartF = StoredExtraEnergyPartF(pot)
      oldextraEnergyPartB = StoredExtraEnergyPartB(pot)

      oldextraEnergy = StoredExtraEnergy(pot) 
      oldextraEnergyF = StoredExtraEnergyF(pot)
      oldextraEnergyB = StoredExtraEnergyB(pot)

      OldTotalEnergy = StoredTotalEnergy(pot) 
      OldTotalEnergyF = StoredTotalEnergyF(pot) 
      OldTotalEnergyB = StoredTotalEnergyB(pot)

      return
      end
