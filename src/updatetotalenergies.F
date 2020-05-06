      subroutine updateTotalEnergies
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'energies.inc'
      include 'solutes.inc'
      include 'flags.inc'
      include 'proteins.inc'
c############################################################
c
c     This routine updates the total energies of the
c     system after a move which involves calls to 
c     'totalEnergy'
c
c     (C) Christopher Woods   29 April 2003
c
c############################################################

      integer i,j

      if (NProteins.gt.0) then
        do i=1,NProteins
          OldProBnd(i)   = NewProBnd(i) 
          OldProAng(i)   = NewProAng(i) 
          OldProDih(i)   = NewProDih(i) 
          OldProUBy(i)   = NewProUBy(i)
          OldProNBC(i)   = NewProNBC(i) 
          OldProNBLJ(i)  = NewProNBLJ(i) 
          OldProBndF(i)  = NewProBndF(i) 
          OldProAngF(i)  = NewProAngF(i) 
          OldProDihF(i)  = NewProDihF(i) 
          OldProUByF(i)  = NewProUByF(i)
          OldProNBCF(i)  = NewProNBCF(i) 
          OldProNBLJF(i) = NewProNBLJF(i)
          OldProBndB(i)  = NewProBndB(i) 
          OldProAngB(i)  = NewProAngB(i) 
          OldProDihB(i)  = NewProDihB(i) 
          OldProUByB(i)  = NewProUByB(i)
          OldProNBCB(i)  = NewProNBCB(i) 
          OldProNBLJB(i) = NewProNBLJB(i)
        enddo

        if (NProteins.gt.1) then
          do i=1,NProteins-1
            do j=i+1,NProteins
              OldProProC(i,j)    = NewProProC(i,j) 
              OldProProLJ(i,j)   = NewProProLJ(i,j)
              OldProProC(j,i)    = NewProProC(j,i) 
              OldProProLJ(j,i)   = NewProProLJ(j,i)
              OldProProCF(i,j)   = NewProProCF(i,j)
              OldProProLJF(i,j)  = NewProProLJF(i,j) 
              OldProProCF(j,i)   = NewProProCF(j,i)
              OldProProLJF(j,i)  = NewProProLJF(j,i) 
              OldProProCB(i,j)   = NewProProCB(i,j)
              OldProProLJB(i,j)  = NewProProLJB(i,j) 
              OldProProCB(j,i)   = NewProProCB(j,i)
              OldProProLJB(j,i)  = NewProProLJB(j,i) 
              OldProProDist(j,i) = NewProProDist(j,i)
              OldProProDist(i,j) = NewProProDist(i,j)
            enddo
          enddo
        endif
      endif
      
      if (NSolutes.gt.0) then
        do i=1,NSolutes
          if (.not.IgnoreSolute(i)) then
            OldSolBnd(i)   = NewSolBnd(i)
            OldSolAng(i)   = NewSolAng(i)
            OldSolDih(i)   = NewSolDih(i)
            OldSolUBy(i)   = NewSolUBy(i)
            OldSolNBC(i)   = NewSolNBC(i)
            OldSolNBLJ(i)  = NewSolNBLJ(i) 
            OldSolBndF(i)  = NewSolBndF(i) 
            OldSolAngF(i)  = NewSolAngF(i) 
            OldSolDihF(i)  = NewSolDihF(i) 
            OldSolUByF(i)  = NewSolUByF(i)
            OldSolNBCF(i)  = NewSolNBCF(i) 
            OldSolNBLJF(i) = NewSolNBLJF(i)
            OldSolBndB(i)  = NewSolBndB(i) 
            OldSolAngB(i)  = NewSolAngB(i) 
            OldSolDihB(i)  = NewSolDihB(i) 
            OldSolUByB(i)  = NewSolUByB(i)
            OldSolNBCB(i)  = NewSolNBCB(i) 
            OldSolNBLJB(i) = NewSolNBLJB(i)
          endif
        enddo

        if (NSolutes.gt.1) then
          do i=1,NSolutes-1
            if (.not.IgnoreSolute(i)) then
              do j=i+1,NSolutes
                if (.not.IgnoreSolute(j)) then 
                  OldSolSolC(j,i)    = NewSolSolC(j,i) 
                  OldSolSolLJ(j,i)   = NewSolSolLJ(j,i)
                  OldSolSolC(i,j)    = NewSolSolC(i,j)  
                  OldSolSolLJ(i,j)   = NewSolSolLJ(i,j)
                  OldSolSolCF(j,i)   = NewSolSolCF(j,i)
                  OldSolSolLJF(j,i)  = NewSolSolLJF(j,i) 
                  OldSolSolCF(i,j)   = NewSolSolCF(i,j)
                  OldSolSolLJF(i,j)  = NewSolSolLJF(i,j) 
                  OldSolSolCB(j,i)   = NewSolSolCB(j,i)
                  OldSolSolLJB(j,i)  = NewSolSolLJB(j,i) 
                  OldSolSolCB(i,j)   = NewSolSolCB(i,j)
                  OldSolSolLJB(i,j)  = NewSolSolLJB(i,j) 
                  OldSolSolDist(j,i) = NewSolSolDist(j,i)
                  OldSolSolDist(i,j) = NewSolSolDist(i,j)
                endif
              enddo
            endif
          enddo
        endif
      
        if (NProteins.gt.0) then
          do i=1,NProteins
            do j=1,NSolutes
              if (.not.IgnoreSolute(j)) then
                OldProSolC(i,j)    = NewProSolC(i,j)           
                OldProSolLJ(i,j)   = NewProSolLJ(i,j)
                OldProSolCF(i,j)   = NewProSolCF(i,j)
                OldProSolLJF(i,j)  = NewProSolLJF(i,j) 
                OldProSolCB(i,j)   = NewProSolCB(i,j)
                OldProSolLJB(i,j)  = NewProSolLJB(i,j) 
                OldProSolDist(i,j) = NewProSolDist(i,j)
              endif
            enddo
          enddo
        endif
      endif
      
      if (NSolvents.gt.0) then
        OldSvnSvnC    = NewSvnSvnC 
        OldSvnSvnLJ   = NewSvnSvnLJ
        OldSvnSvnCF   = NewSvnSvnCF
        OldSvnSvnLJF  = NewSvnSvnLJF 
        OldSvnSvnCB   = NewSvnSvnCB
        OldSvnSvnLJB  = NewSvnSvnLJB 
        OldSvnSvnDist = NewSvnSvnDist
      
        if (NProteins.gt.0) then
          do i=1,NProteins
            OldProSvnC(i)    = NewProSvnC(i) 
            OldProSvnLJ(i)   = NewProSvnLJ(i)
            OldProSvnCF(i)   = NewProSvnCF(i)
            OldProSvnLJF(i)  = NewProSvnLJF(i) 
            OldProSvnCB(i)   = NewProSvnCB(i)
            OldProSvnLJB(i)  = NewProSvnLJB(i) 
            OldProSvnDist(i) = NewProSvnDist(i)
          enddo
        endif
        
        if (NGCSolutes.gt.0) then
          do i=1,NGCSolutes
              OldIndGCSSvnC(i)    = NewIndGCSSvnC(i)
              OldIndGCSSvnLJ(i)   = NewIndGCSSvnLJ(i)
              OldIndGCSSvnCF(i)   = NewIndGCSSvnCF(i)
              OldIndGCSSvnLJF(i)  = NewIndGCSSvnLJF(i)
              OldIndGCSSvnCB(i)   = NewIndGCSSvnCB(i)
              OldIndGCSSvnLJB(i)  = NewIndGCSSvnLJB(i)
              OldGCSSvnDist(i)    = NewGCSSvnDist(i)
         enddo
        endif
      
        if (NSolutes.gt.0) then
          do i=1,NSolutes
            if (.not.IgnoreSolute(i)) then
              OldSolSvnC(i)    = NewSolSvnC(i)  
              OldSolSvnLJ(i)   = NewSolSvnLJ(i)
              OldSolSvnCF(i)   = NewSolSvnCF(i)
              OldSolSvnLJF(i)  = NewSolSvnLJF(i) 
              OldSolSvnCB(i)   = NewSolSvnCB(i)
              OldSolSvnLJB(i)  = NewSolSvnLJB(i) 
              OldSolSvnDist(i) = NewSolSvnDist(i)
              if (i.eq.SvnPrefSol .and. PrefSample) then
c               save the solvent probabilities...
                do j=1,NSolvents
                  SvnProb(j) = SvnNewProb(j)
                enddo
                SvnSumProb = SvnNewSumProb
                SvnMaxProb = SvnNewMaxProb
              endif
            endif
          enddo
        endif
      endif
     
      if (NGCSolutes.eq.0) then
        OldIndGCSC = NewIndGCSC
        OldIndGCSLJ = NewIndGCSLJ
        OldIndGCSCF = NewIndGCSCF
        OldIndGCSLJF = NewIndGCSLJF
        OldIndGCSCB = NewIndGCSCB
        OldIndGCSLJB = NewIndGCSLJB
        OldGCSGCSC  = NewGCSGCSC
        OldGCSGCSLJ  = NewGCSGCSLJ
        OldGCSGCSCF  = NewGCSGCSCF
        OldGCSGCSLJF  = NewGCSGCSLJF
        OldGCSGCSCB  = NewGCSGCSCB
        OldGCSGCSLJB  = NewGCSGCSLJB

        if (NProteins.gt.0) then
          do i=1,NProteins
            OldIndGCSProC(i)    = NewIndGCSProC(i)
            OldIndGCSProLJ(i)   = NewIndGCSProLJ(i)
            OldIndGCSProCF(i)   = NewIndGCSProCF(i)
            OldIndGCSProLJF(i)  = NewIndGCSProLJF(i)
            OldIndGCSProCB(i)   = NewIndGCSProCB(i)
            OldIndGCSProLJB(i)  = NewIndGCSProLJB(i)
          enddo
        endif
      endif
      
      OldTotalEnergy  = NewTotalEnergy 
      OldTotalEnergyF = NewTotalEnergyF
      OldTotalEnergyB = NewTotalEnergyB
            
      return
      end
