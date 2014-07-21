c##################################
c Computes the Generalised Born energy 
c of the *entire* system
c     Julien Michel 2004
c##################################

      double precision function GBEnergy(savenrg, gbnrgf, gbnrgb)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'printer.inc'
      include 'simparams.inc'
      include 'enums.inc'
      include 'flags.inc'
      include 'gb.inc'
      include 'quickpot.inc'
      integer i,j,k,l,savenrg,nres,nresk
      double precision nrg,gbnrgf,gbnrgb
      double precision gbnrgintraf,gbnrgintrab,gbnrginterf,gbnrginterb
      double precision ResidueIntraGBEnergy
      double precision ResidueResidueInterGBEnergy
      double precision SoluteIntraGBEnergy
      double precision SoluteSoluteInterGBEnergy
      double precision ResidueSoluteGBEnergy
      double precision SolventIntraGBEnergy
      double precision SolventSolventInterGBEnergy
      double precision ResidueSolventGBEnergy
      double precision SoluteSolventGBEnergy
c      integer nats,resstart,resend
      double precision scl,sclf,sclb
      double precision getLowScale,getHighScale

c      write(printstring,*) 'GB ENERGY WAS CALLED ',savenrg
c      call printLine(DEBUG,printstring)

c     Prepare the calculation
      GBEnergy = OldGBEnergy
      gbnrgf = OldGBEnergyF
      gbnrgb = OldGBEnergyB
      nrg = ZERO
c     Bad things happens if savenrg is not OLD or NEW !
      if (savenrg .eq. OLD) then
        if (.not. ConstBorn) then 
          call ZeroPDACont
          call TotalBorn(savenrg)
          call pushBackPDACont
        endif
c       Save coordinates of the parts that have moved, this is so that
c       when we want to calculate the gb energy of a new configuration
c       the subroutines 'markgbupdate' knows which residues have moved
        call gbsavesystem        
c       Return as we now what the energy is 
        return
      else if (savenrg .eq. NEW ) then
        if (.not. ConstBorn) then 
          call TotalBorn(savenrg)
          call GetEffectiveBornRadii
        endif
c       Flag Moving bits
        call markgbupdate        
      endif

c      write(printstring,*) 'STARTING WITH OLDGB ',GBEnergy,' gbnrgf ',gbnrgf,
c     &  ' gbnrgb ',gbnrgb
c      call printLine(DEBUG,printstring)

c
c Now do the GB Energy
c      
      if (NProteins.gt.0) then      
        do 100 i=1,NProteins
          nres = NResidues(i)
          do 101 j=1,nres
c         Do pair pair intra protein    
c            write (printstring,*) 'DOING FAST PROT ',i,' RES INTRA ',j
c            call printLine(DEBUG,printstring)  
            if (calcResidueGB(i,j)) then
              nrg = ResidueIntraGBEnergy(i,j,gbnrgintraf,gbnrgintrab)
c              write(printstring,*) 'PROT ',i,' RES INTRA ',j,
c     &          ' NRG CHANGE ',nrg, ' Fwd ',gbnrgintraf
c              call printLine(DEBUG,printstring)
              GBEnergy = GBEnergy + nrg
              gbnrgf = gbnrgf + gbnrgintraf
              gbnrgb = gbnrgb + gbnrgintrab
            endif
c           if (savenrg .eq. OLD) then
c             OldProGBIntra(i,j) = nrg
c             OldProGBIntraF(i,j) = gbnrgintraf
c             OldProGBIntraB(i,j) = gbnrgintrab
c           else if (savenrg .eq. NEW) then
c             NewProGBIntra(i,j) = nrg
c             NewProGBIntraF(i,j) = gbnrgintraf
c             NewProGBIntraB(i,j) = gbnrgintrab
c           endif
 101      continue

c        write(printstring,*) 'NRG AFTER PROT INTRA',FastGBEnergy
c        call printLine(DEBUG,printstring)
c       Now the contributions between residues in the same protein
          if (nres .gt. 1) then
            do 102 j=1,nres-1
              do 103 k=j+1,nres
c                write (printstring,*) 'DOING FAST PROT ',i,' RES ',j,
c     &            ' RES ',k
c                call printLine(DEBUG,printstring)   
                if (calcResidueGB(i,j) .or. calcResidueGB(i,k) ) then
                  nrg = ResidueResidueInterGBEnergy(i,j,i,k,gbnrginterf,gbnrginterb)
c                write(printstring,*) 'PROT ',i,' INTER RES ',j,' RES ',k,
c     &            'NRG CHANGE ',nrg
c                call printLine(DEBUG,printstring)
                GBEnergy = GBEnergy + nrg
                gbnrgf = gbnrgf + gbnrginterf
                gbnrgb = gbnrgb + gbnrginterb
              endif
 103          continue
 102        continue
          endif
 100    continue
c       write(printstring,*) 'NRG AFTER PROT INTRA INTER',FastGBEnergy
c       call printLine(DEBUG,printstring)
        if (NProteins .gt. 1) then
c         Do the contributions between a residue and all the residues in every 
c         other protein
          do 104 i=1,NProteins-1
            nres = NResidues(i)
            do 105 j=1,nres
              do 106 k=i+1,NProteins
                nresk = NResidues(k)
                do 107 l=1,nresk
c                 write (printstring,*) 'DOING FAST PROT ',i,' RES ',j,
c     &            ' PROT ',k,' RES ',l
c                 call printLine(DEBUG,printstring) 
                  if (calcResidueGB(i,j) .or. calcResidueGB(k,l)) then
                    nrg = ResidueResidueInterGBEnergy(i,j,k,l,gbnrginterf,gbnrginterb)
                    GBEnergy = GBEnergy + nrg
                    gbnrgf = gbnrgf + gbnrginterf
                    gbnrgb = gbnrgb + gbnrginterb
                  endif
 107            continue
 106          continue
 105        continue
 104      continue
        endif
      endif
c      write(printstring,*) 'NRG AFTER PROT PROT INTER ',FastGBEnergy
c      call printLine(DEBUG,printstring)
      if (NSolutes.gt.0) then
        do 200 i=1,NSolutes
          if (IgnoreSolute(i)) goto 200
c         Do pair-pair intrasolute
          if (calcSoluteGB(i)) then
            nrg = SoluteIntraGBEnergy(i,gbnrgintraf,gbnrgintrab) 
c           see if the solute is part of a dual topology solute
            if (DualTopologySolute(i).gt.0) then
c             this is the low dual topology solute
              scl = getLowScale(Lambda)
              sclb = getLowScale(LambdaB)
              sclf = getLowScale(LambdaF)
            else if (DualTopologySolute(i).eq.0) then
c             this is not a dual topology solute
              scl = ONE
              sclf = ONE
              sclb = ONE
            else
c             this is the high dual topology solute
              scl = getHighScale(Lambda)
              sclb = getHighScale(LambdaB)
              sclf = getHighScale(LambdaF)
            endif
c            write(printstring,*) 'Solute ',i,' Change GBintra ',nrg*scl,
c     &        ' F ',gbnrgintraf*sclf,' B ',gbnrgintrab*sclb
c            call printLine(DEBUG,printstring)
            GBEnergy = GBEnergy + nrg * scl
            gbnrgf = gbnrgf + gbnrgintraf * sclf
            gbnrgb = gbnrgb + gbnrgintrab * sclb
c            if (savenrg .eq. OLD) then
c              OldSolGBIntra(i) = nrg
c              OldSolGBIntraF(i) = gbnrgintraf
c              OldSolGBIntraB(i) = gbnrgintrab
c            else if (savenrg .eq. NEW) then
c              NewSolGBIntra(i) = nrg
c              NewSolGBIntraF(i) = gbnrgintraf
c              NewSolGBIntraB(i) = gbnrgintrab
c            endif
          endif
 200    continue
c        write(printstring,*) 'NRG AFTER SOL INTRA',FastGBEnergy
c        call printLine(DEBUG,printstring)
        if (NSolutes.gt.1) then
          do 201 i=1,NSolutes-1
            if (IgnoreSolute(i)) goto 201
            do 202 j=i+1,NSolutes
              if (IgnoreSolute(j)) goto 202
c             Do pair-pair inter solutes
c              write (printstring,*) 'DOING SOLUTE SOLUTE ',i,'-',j
c             call printLine(DEBUG,printstring)
              if (calcSoluteGB(i) .or. calcSoluteGB(j) ) then
                nrg = SoluteSoluteInterGBEnergy(i,j,gbnrginterf,gbnrginterb)
c                write(printstring,*) 'i ',i,' j ',j,
c     &            ' EGBinter is ',nrg,' F ',gbnrginterf,' B ',gbnrginterb
c                call printLine(DEBUG,printstring)
                GBEnergy = GBEnergy + nrg
                gbnrgf = gbnrgf + gbnrginterf
                gbnrgb = gbnrgb + gbnrginterb
              endif
 202        continue
 201      continue
        endif
c       write(printstring,*) 'NRG AFTER SOL-SOL INTER',FastGBEnergy
c       call printLine(DEBUG,printstring)
c       
c       Now do the solutes proteins contributions
c       
        if (NProteins.gt.0) then
          do 210 i=1,NProteins
            nres = NResidues(i)
            do 220 j=1,nres
              do 230 k=1,NSolutes
                if (IgnoreSolute(k)) goto 230
c                write (printstring,*) 'DOING FAST PROT ',i,' RES ',j,
c     &            ' SOLUTE ',k
c                call printLine(DEBUG,printstring)
                if (calcResidueGB(i,j) .or. calcSoluteGB(k)) then
                  nrg = ResidueSoluteGBEnergy(i,j,k,gbnrginterf,gbnrginterb)
c                 see if this is a dual topology solute
                  if (DualTopologySolute(k).gt.0) then
                    scl = getLowScale(Lambda)
                    sclb = getLowScale(LambdaB)
                    sclf = getLowScale(LambdaF)
                  else if (DualTopologySolute(k).eq.0) then
                    scl = ONE
                    sclb = ONE
                    sclf = ONE
                  else
                    scl = getHighScale(Lambda)
                    sclb = getHighScale(LambdaB)
                    sclf = getHighScale(LambdaF)
                  endif
c                 write(printstring,*) 'NRG CHANGE ',nrg
c                 call printLine(DEBUG,printstring)
                  if (SoluteSoftCore(k)) then 
                    nrg = nrg * scl**GBSoftPower
                    gbnrginterf = gbnrginterf * sclf**GBSoftPower
                    gbnrginterb = gbnrginterb * sclb**GBSoftPower
                  endif
                  GBEnergy = GBEnergy + nrg * scl
                  gbnrgf = gbnrgf + gbnrginterf * sclf
                  gbnrgb = gbnrgb + gbnrginterb * sclb
                endif
 230          continue
 220        continue
 210      continue
        endif
      endif

c       now do the same for the solvent
      if (NSolvents.gt.0) then
        do 300 i=1,NSolvents
c         write (printstring,*) 'DOING SOLVENT INTRA ',i
c         call printLine(DEBUG,printstring)
          if (calcSolventGB(i)) then
            nrg = SolventIntraGBEnergy(i,gbnrgintraf,gbnrgintrab)
c           write (printstring,*) i,' NRG CHANGE SOLV INTRA ',nrg
c           call printLine(DEBUG,printstring)
            GBEnergy = GBEnergy + nrg 
            gbnrgf = gbnrgf + gbnrgintraf 
            gbnrgb = gbnrgb + gbnrgintrab
          endif
 300    continue

        if (NSolvents .gt. 1) then
          do 301 i=1,NSolvents-1
            do 302 j=i+1,NSolvents
c             write (printstring,*) 'DOING SOLVENT SOLVENT ',i,'-',j
c             call printLine(DEBUG,printstring)
              if (calcSolventGB(i) .or. calcSolventGB(j)) then
                nrg = SolventSolventInterGBEnergy(i,j,gbnrginterf,gbnrginterb)
c               write (printstring,*) i,'-',j,' NRG CHANGE SOLV INTER ',nrg
c               call printLine(DEBUG,printstring)
                GBEnergy = GBEnergy + nrg
                gbnrgf = gbnrgf + gbnrginterf
                gbnrgb = gbnrgb + gbnrginterb
              endif
 302        continue
 301      continue
        endif
c       Now the solvent proteins 
c        
        if (NProteins.gt.0) then
          do 310 i=1,NProteins
            nres = NResidues(i)
            do 311 j=1,nres
              do 312 k=1,NSolvents
c               write (printstring,*) 'DOING PROT ',i,' RES ',j,
c               &            ' SOLVENT ',k
c               call printLine(DEBUG,printstring)
                if (calcResidueGB(i,j) .or. calcSolventGB(k)) then
                  nrg = ResidueSolventGBEnergy(i,j,k,gbnrginterf,gbnrginterb)
                  GBEnergy = GBEnergy + nrg
                  gbnrgf = gbnrgf + gbnrginterf
                  gbnrgb = gbnrgb + gbnrginterb
                endif
 312          continue
 311        continue
 310      continue
        endif
c       Do the solute solvent contributions
        if (NSolutes .gt. 0) then
          do 320 i=1,NSolutes
            if (IgnoreSolute(i)) goto 320
c           see if this is a dual topology solute
            if (DualTopologySolute(i).gt.0) then
              scl = getLowScale(Lambda)
              sclb = getLowScale(LambdaB)
              sclf = getLowScale(LambdaF)
            else if (DualTopologySolute(i).eq.0) then
              scl = ONE
              sclb = ONE
              sclf = ONE
            else
              scl = getHighScale(Lambda)
              sclb = getHighScale(LambdaB)
              sclf = getHighScale(LambdaF)
            endif
            do 330 j=1,NSolvents
c             write (printstring,*) 'DOING SOLUTE ',i,' SOLVENT ',j
c             call printLine(DEBUG,printstring)
              if (calcSoluteGB(i) .or. calcSolventGB(j)) then
                nrg = SoluteSolventGBEnergy(i,j,gbnrginterf,gbnrginterb)
                if (SoluteSoftCore(i)) then
                  nrg = nrg * scl**GBSoftPower
                  gbnrginterf = gbnrginterf * sclf**GBSoftPower
                  gbnrginterb = gbnrginterb * sclb**GBSoftPower
                endif
                GBEnergy = GBEnergy + nrg * scl
                gbnrgf = gbnrgf + gbnrginterf * sclf
                gbnrgb = gbnrgb + gbnrginterb * sclb
              endif
 330        continue
 320      continue
        endif
      endif
      
      if (savenrg .eq. OLD) then
c       We never should get here
        OldGBEnergy = GBEnergy
        OldGBEnergyF = gbnrgf
        OldGBEnergyB = gbnrgb
      else if (savenrg .eq. NEW) then
        NewGBEnergy = GBEnergy
        NewGBEnergyF = gbnrgf
        NewGBEnergyB = gbnrgb
      endif

c      write(printstring,*) 'FAST GB ',GBEnergy,' Fwd ',gbnrgf,
c     &  ' Bwd ',gbnrgb
c      call printLine(DEBUG,printstring)

c
c     Clear bits saved on the lists
c
      call clearSaveSoluteGB
      call clearSaveSolventGB
      call clearSaveResidueGB

      return 
      end
