      subroutine StoreEnergies(pot)
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
c############################################################
c     This subroutine save all the energy components 
c     of the configuration in the 'OLD' state
c
c
c     (c) Julien Michel April 2005
c###########################################################
      integer i,j,k,pot

c      write(printstring,*) 'STORING ENERGIES FOR POT ',pot
c      call printLine(DEBUG,printstring)

      if (NProteins .gt. 0) then
        do i=1,NProteins
          StoredProBnd(i,pot) = OldProBnd(i)
          StoredProAng(i,pot) = OldProAng(i)
          StoredProDih(i,pot) = OldProDih(i) 
          StoredProUBy(i,pot) = OldProUBy(i) 
          StoredProNBC(i,pot) = OldProNBC(i) 
          StoredProNBLJ(i,pot) = OldProNBLJ(i) 
          StoredProBndF(i,pot) = OldProBndF(i)
          StoredProAngF(i,pot) = OldProAngF(i)
          StoredProDihF(i,pot) = OldProDihF(i)
          StoredProUByF(i,pot) = OldProUByF(i) 
          StoredProNBCF(i,pot) = OldProNBCF(i) 
          StoredProNBLJF(i,pot) = OldProNBLJF(i)
          StoredProBndB(i,pot) = OldProBndB(i) 
          StoredProAngB(i,pot) = OldProAngB(i) 
          StoredProDihB(i,pot) = OldProDihB(i)
          StoredProUByB(i,pot) = OldProUByB(i)
          StoredProNBCB(i,pot) = OldProNBCB(i)
          StoredProNBLJB(i,pot) = OldProNBLJB(i)
        enddo
        if (NProteins.gt.1) then
          do i=1,NProteins-1
            do j=i+1,NProteins
              StoredProProC(i,j,pot) = OldProProC(i,j)
              StoredProProLJ(i,j,pot) = OldProProLJ(i,j) 
              StoredProProC(j,i,pot) = OldProProC(j,i) 
              StoredProProLJ(j,i,pot) = OldProProLJ(j,i) 
              StoredProProCF(i,j,pot) = OldProProCF(i,j) 
              StoredProProLJF(i,j,pot) = OldProProLJF(i,j) 
              StoredProProCF(j,i,pot) = OldProProCF(j,i) 
              StoredProProLJF(j,i,pot) = OldProProLJF(j,i)
              StoredProProCB(i,j,pot) = OldProProCB(i,j)
              StoredProProLJB(i,j,pot) = OldProProLJB(i,j)
              StoredProProCB(j,i,pot) = OldProProCB(j,i)
              StoredProProLJB(j,i,pot) = OldProProLJB(j,i) 
              StoredProProDist(j,i,pot) = OldProProDist(j,i) 
              StoredProProDist(i,j,pot) = OldProProDist(i,j) 
            enddo
          enddo
        endif
      endif

      if (NSolutes.gt.0) then
        do i=1,NSolutes
          if (.not.IgnoreSolute(i)) then
            StoredSolBnd(i,pot) = OldSolBnd(i)
            StoredSolAng(i,pot) = OldSolAng(i)
            StoredSolDih(i,pot) = OldSolDih(i) 
            StoredSolUBy(i,pot) = OldSolUBy(i)
            StoredSolNBC(i,pot) = OldSolNBC(i)
            StoredSolNBLJ(i,pot) = OldSolNBLJ(i)
            StoredSolBndF(i,pot) = OldSolBndF(i)
            StoredSolAngF(i,pot) = OldSolAngF(i)
            StoredSolDihF(i,pot) = OldSolDihF(i)
            StoredSolUByF(i,pot) = OldSolUByF(i)
            StoredSolNBCF(i,pot) = OldSolNBCF(i)
            StoredSolNBLJ(i,pot) = OldSolNBLJF(i)
            StoredSolBndB(i,pot) = OldSolBndB(i)
            StoredSolAngB(i,pot) = OldSolAngB(i)
            StoredSolDihB(i,pot) = OldSolDihB(i)
            StoredSolUByB(i,pot) = OldSolUByB(i)
            StoredSolNBCB(i,pot) = OldSolNBCB(i)
            StoredSolNBLJB(i,pot) = OldSolNBLJB(i)
          endif
        enddo
        
        if (NSolutes.gt.1) then
          do i=1,NSolutes-1
            if (.not.IgnoreSolute(i)) then
              do j=i+1,NSolutes
                if (.not.IgnoreSolute(j)) then 
                  StoredSolSolC(j,i,pot) = OldSolSolC(j,i)
                  StoredSolSolLJ(j,i,pot) = OldSolSolLJ(j,i)
                  StoredSolSolC(i,j,pot) = OldSolSolC(i,j)
                  StoredSolSolLJ(i,j,pot) = OldSolSolLJ(i,j)
                  StoredSolSolCF(j,i,pot) = OldSolSolCF(j,i) 
                  StoredSolSolLJF(j,i,pot) = OldSolSolLJF(j,i)
                  StoredSolSolCF(i,j,pot) = OldSolSolCF(i,j) 
                  StoredSolSolLJF(i,j,pot) = OldSolSolLJF(i,j) 
                  StoredSolSolCB(j,i,pot) = OldSolSolCB(j,i) 
                  StoredSolSolLJB(j,i,pot) = OldSolSolLJB(j,i)
                  StoredSolSolCB(i,j,pot) = OldSolSolCB(i,j) 
                  StoredSolSolLJB(i,j,pot) = OldSolSolLJB(i,j) 
                  StoredSolSolDist(j,i,pot) = OldSolSolDist(j,i)
                  StoredSolSolDist(i,j,pot) = OldSolSolDist(i,j) 
                endif
              enddo
            endif
          enddo
        endif

        if (NProteins.gt.0) then
          do j=1,NSolutes   
            if (.not.IgnoreSolute(j)) then
              do i=1,NProteins
                StoredProSolC(i,j,pot) = OldProSolC(i,j)
                StoredProSolLJ(i,j,pot) = OldProSolLJ(i,j)
                StoredProSolCF(i,j,pot) = OldProSolCF(i,j) 
                StoredProSolLJF(i,j,pot) = OldProSolLJF(i,j)
                StoredProSolCB(i,j,pot) = OldProSolCB(i,j) 
                StoredProSolLJB(i,j,pot) = OldProSolLJB(i,j)
                StoredProSolDist(i,j,pot) = OldProSolDist(i,j)
              enddo
            endif
          enddo
        endif

      endif

      if (NSolvents .gt. 0) then
        StoredSvnSvnC(pot) = OldSvnSvnC
        StoredSvnSvnLJ(pot) = OldSvnSvnLJ
        StoredSvnSvnCF(pot) = OldSvnSvnCF
        StoredSvnSvnLJF(pot) = OldSvnSvnLJF
        StoredSvnSvnCB(pot) = OldSvnSvnCB
        StoredSvnSvnLJB(pot) = OldSvnSvnLJB
        StoredSvnSvnDist(pot) = OldSvnSvnDist
        if (NProteins.gt.0) then
          do i=1,NProteins
            StoredProSvnC(i,pot) = OldProSvnC(i) 
            StoredProSvnLJ(i,pot) = OldProSvnLJ(i) 
            StoredProSvnCF(i,pot) = OldProSvnCF(i) 
            StoredProSvnLJF(i,pot) = OldProSvnLJF(i)
            StoredProSvnCB(i,pot) = OldProSvnCB(i) 
            StoredProSvnLJB(i,pot) = OldProSvnLJB(i)
            StoredProSvnDist(i,pot) = OldProSvnDist(i)
          enddo
        endif
        if (NSolutes.gt.0) then
          do i=1,NSolutes
            if (.not.IgnoreSolute(i)) then
              StoredSolSvnC(i,pot) = OldSolSvnC(i) 
              StoredSolSvnLJ(i,pot) = OldSolSvnLJ(i) 
              StoredSolSvnCF(i,pot) = OldSolSvnCF(i) 
              StoredSolSvnLJF(i,pot) = OldSolSvnLJF(i) 
              StoredSolSvnCB(i,pot) = OldSolSvnCB(i) 
              StoredSolSvnLJB(i,pot) = OldSolSvnLJB(i) 
              StoredSolSvnDist(i,pot) = OldSolSvnDist(i) 
              if (i.eq.SvnPrefSol .and. PrefSample) then
c     save the solvent probabilities...
                 do j=1,NSolvents
                    StoredSvnProb(j,pot) = SvnProb(j) 
                 enddo
                 StoredSumProb(pot) = SvnSumProb 
                 StoredMaxProb(pot) = SvnMaxProb
              endif
c     save the residue probabilities
              if (i.eq.ProPrefSol .and. ProPrefSample) then
                 do j=1,NProteins
                    do k=1,NResidues(j)
                       StoredProProb(j,k,pot) = ProProb(j,k)
                    enddo
                 enddo
                 StoredProSumProb(pot) = ProSumProb
                 StoredProMaxProb(pot) = ProMaxProb
              endif
              
            endif
          enddo
        endif
      endif
c
c     Need to store the extra energies
c     
c     
      if (SolventCap) then 
        StoredSolventCapPart(pot) = oldSvnCapEnergyPart
        StoredSolventCap(pot) = oldSvnCapEnergy
      endif

      if (NHarmonicRestraints .gt. 0) then 
        StoredHarmEnergyPart(pot) = oldHarmEnergyPart
        StoredHarmEnergyPartB(pot) = oldHarmEnergyPartB
        StoredHarmEnergyPartF(pot) = oldHarmEnergyPartF
        StoredHarmEnergy(pot) = oldHarmEnergy
        StoredHarmEnergyB(pot) = oldHarmEnergyB
        StoredHarmEnergyF(pot) = oldHarmEnergyF
      endif

      if (NFlatBottomRestraints .gt. 0) then
        StoredFBotEnergyPart(pot) = oldFBotEnergyPart
        StoredFBotEnergyPartB(pot) = oldFBotEnergyPartB
        StoredFBotEnergyPartF(pot) = oldFBotEnergyPartF
        StoredFBotEnergy(pot) = oldFBotEnergy
        StoredFBotEnergyB(pot) = oldFBotEnergyB
        StoredFBotEnergyF(pot) = oldFBotEnergyF
      endif

      if (doSA) then
        StoredTotalSurface(pot) = oldTotalSurface
        StoredTotalSurfaceF(pot) = oldTotalSurfaceF 
        StoredTotalSurfaceB(pot) = oldTotalSurfaceB 
        StoredSAEnergy(pot) = oldSAEnergy 
        StoredSAEnergyF(pot) = oldSAEnergyF
        StoredSAEnergyB(pot) = oldSAEnergyB 
      endif

      if (doGB) then
c        write(printstring,*) 'STORING GB FOR POT ',pot,' ',OldGBEnergy
c        call printLine(DEBUG,printstring)
        StoredGBEnergy(pot) = OldGBEnergy 
        StoredGBEnergyF(pot) = OldGBEnergyF 
        StoredGBEnergyB(pot) = OldGBEnergyB
      endif

      StoredExtraEnergyPart(pot) = oldextraEnergyPart
      StoredExtraEnergyPartF(pot) = oldextraEnergyPartF
      StoredExtraEnergyPartB(pot) = oldextraEnergyPartB

      StoredExtraEnergy(pot) = oldextraEnergy
      StoredExtraEnergyF(pot) = oldextraEnergyF
      StoredExtraEnergyB(pot) = oldextraEnergyB

      StoredTotalEnergy(pot) = OldTotalEnergy
      StoredTotalEnergyF(pot) = OldTotalEnergyF
      StoredTotalEnergyB(pot) = OldTotalEnergyB 

      return
      end
