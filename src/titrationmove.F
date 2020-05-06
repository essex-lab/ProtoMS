      
      logical function titrationMove(updatenrg,movnum,movtyp)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'flags.inc'
      include 'moves.inc'
      include 'enums.inc'
      include 'energies.inc'
c#################################################################
c
c     This function performs everything to do with a titration
c     move... - returns true if the move was accepted
c
c     If updatenrg = 0 then don't bother updating the energy
c     components during the move. Else do keep updating the
c     energy components. (use 0 for equilibration, 1 for collection)
c
c     (C) Christopher Woods  4 February 2002
c
c#################################################################

      integer iprot1,iprot2,ires1,ires2,movnum,movtyp,updatenrg
      integer iacid,ibase,oldprot,savenrg
      integer i
      logical protonate,deprotonate
      logical metropolis
      logical convertRes
      logical flipHistidine
      integer ihis(MAXBASICRES),nhis
      double precision ieold,ienew,nrgb,nrgf
      double precision residueEnergy
      double precision extraEnergy,inrgf,inrgb
      character*4 statstring
      character*4 oldres1,oldres2
      integer ranInt      
      
c     This subroutine randomly chooses two complementary
c     tiratable groups on the protein, and moves a proton
c     between them (eventually, this routine will use the
c     'pH' of the solvent to add and take protons as well!)
      titrationMove = .false.
      if (NProteins.le.0) return

c     clear the move lists
      call clearMoveList()      
            
c     make sure that updatenrg is valid
      if (updatenrg.eq.UPDATEENERGY) then
        savenrg = OLD
      else
        savenrg = DONTCARE
      endif

c     the first group - randomly choose the protein...
      if (NProteins.eq.1) then
        iprot1 = 1
      else
        iprot1 = ranInt(1,NProteins)
      endif
      
c     there is a 1 in 3 chance that we will only attempt to flip the hydrogen on
c     a histidine
      if (ranInt(1,3).eq.1) then
c       flipping a histidine - first find a histidine residue in iprot1

c       there are no basic residues, and thus no available histidines! - thus continue with a normal
        if (NBasicRes(iprot1).lt.1) goto 99
c       run through all basic residues to create a list of histidines...   
        nhis = 0
        do i=1,NBasicRes(iprot1)
          if (ResNam(iprot1,BasicRes(iprot1,i)).eq.'HIS' .or. ResNam(iprot1,BasicRes(iprot1,i)).eq.'HID') then
            nhis = nhis + 1
            ihis(nhis) = BasicRes(iprot1,i)
          endif
        enddo

c       if there are no histidines, then continue with a normal titration move
        if (nhis.lt.1) goto 99
        
c       now randomly choose one histidine from the list
        ires1 = ihis(ranInt(1,nhis))
        
c       set this residue as having moved
        call setResidueMoved(iprot1,ires1)
        
c       calculate the old energy of the histidine
        ieold = residueEnergy(savenrg,iprot1,ires1,nrgb,nrgf) + 
     .             extraEnergy(savenrg,inrgb,inrgf)
        oldres1 = ResNam(iprot1,ires1)
c       flip the hydrogen
        if (.not.flipHistidine(iprot1,ires1)) then
          write(printstring,*) "Problem flipping the histidine of ",iprot1," ",ires1
          call closeProgram(-1,printstring)
        endif

c       switch to new energy if necessary
        if (updatenrg.eq.UPDATEENERGY) savenrg = NEW
        
c       calculate the new energy of the histidine
        ienew = residueEnergy(savenrg,iprot1,ires1,nrgb,nrgf) + 
     .            extraEnergy(savenrg,inrgb,inrgf)
        
c       perform a metropolis test...
        if (metropolis(ieold,ienew,movtyp)) then
c         if pass, then leave as it is
          titrationMove = .true.
          statstring = "PASS"
          if (updatenrg.eq.UPDATEENERGY) then
c           update the energies from residue ires2
            call updateResidueEnergies(iprot1)
            call updateExtraEnergies
          endif
          call clearResidueEnergyStack
        else
c         if fail, then retun the residues to thier previous state...
          titrationMove = .false.
          statstring = "FAIL"
          if (.not.flipHistidine(iprot1,ires1)) then
            write(printstring,*) "Problem flipping the histidine of ",iprot1," ",ires1
            call closeProgram(-1,printstring)
          endif
c         clear the residue energy buffer
          call clearResidueEnergyStack
        endif

        call clearMoveList()
        
c       print the status of this move
        write(printstring,8) movnum,iprot1,oldres1,ResNum(iprot1,ires1),
     .              OldTotalEnergy,OldTotalEnergyF,OldTotalEnergyB,statstring
8       format(i8,": TITRATION. HIS-FLIP. Attempt",i2," ",a,"(",i3,") -> ",
     .              f12.3," | ",f12.3,f12.3," ",a4)
        call printLine(MOVE,printstring)
        
c       since the move is complete, return from the function
        return
        
      endif
        
99    continue
      
c     randomly choose an acidic residue in protein iprot1
      if (NAcidicRes(iprot1).le.0) then
        oldprot = iprot1
100     continue
c       need to choose another protein
        iprot1 = iprot1 + 1
        if (iprot1.gt.NProteins) iprot1 = 1
c       if there are no acidic residues then return!
        if (iprot1.eq.oldprot) then
          write(printstring,*) "There are no acidic residues to change!"
          call printLine(WARNING,printstring)
          return
        endif
        if (NAcidicRes(iprot1).le.0) goto 100
      endif
      
c     use a function to generate a random integer between
c     1 and NAcidicRes
      iacid = ranInt(1,NAcidicRes(iprot1))
      ires1 = AcidicRes(iprot1,iacid)
      
c     now randomly choose another protein...
      if (NProteins.eq.1) then
        iprot2 = 1
      else
c       use a function to generate a random integer between
c       1 and NProteins
        iprot2 = ranInt(1,NProteins)
      endif
      
      if (NBasicRes(iprot2).le.0) then
        oldprot = iprot2
20      continue
        iprot2 = iprot2 + 1
        if (iprot2.gt.NProteins) iprot2 = 1
        if (iprot2.eq.oldprot) then 
          write(printstring,*) "There are no basic residues to change!"
          call printLine(WARNING,printstring)
          return
        endif        
        if (NBasicRes(iprot2).le.0) return
      endif
      
c     use a function to calculate a random integer between
c     1 and NBasicRes
      ibase = ranInt(1,NBasicRes(iprot2))
      ires2 = BasicRes(iprot2,ibase)

c     set these two residues as having moved
      call setResidueMoved(iprot1,ires1)
      call setResidueMoved(iprot2,ires2)     
            
c     calculate the energies of the residues as they are...
c     since we are moving two residues, we need to use the history
c     to save the energies from the first residues move
      ieold = residueEnergy(savenrg,iprot1,ires1,nrgb,nrgf) + 
     .          extraEnergy(savenrg,inrgb,inrgf)
     
      DoneResNrg(iprot1,ires1) = .true.
c     save this residue's energy
      call pushBackOldResidueEnergy
      ieold = ieold + residueEnergy(savenrg,iprot2,ires2,nrgb,nrgf)      
      DoneResNrg(iprot2,ires2) = .true.

      DoneResNrg(iprot1,ires1) = .false.
      DoneResNrg(iprot2,ires2) = .false.

c     save the old residue names - allows us to see what happened!
      oldres1 = ResNam(iprot1,ires1)
      oldres2 = ResNam(iprot2,ires2)

c     protonate/deprotonate the residues..
      if (deprotonate(iprot1,ires1)) then
        if (.not.protonate(iprot2,ires2)) then
          write(printstring,10) ResNam(iprot2,ires2),ResNum(iprot2,ires2)
10        format("Problem protonating ",a,"(",i3,")")
          call closeProgram(-1,printstring)
        endif
      else
        write(printstring,11) ResNam(iprot1,ires1),ResNum(iprot1,ires1)
11      format("Problem deprotonating ",a,"(",i3,")")
        call closeProgram(-1,printstring)
      endif
     
c     switch to new energy if necessary
      if (updatenrg.eq.UPDATEENERGY) savenrg = NEW
      
c     calculate the energies of the residues now...
      ienew = residueEnergy(savenrg,iprot1,ires1,nrgb,nrgf)
      DoneResNrg(iprot1,ires1) = .true.
c     save ires1's new residue energy
      call pushBackNewResidueEnergy
c     now calculate the new energy of residue 2
      ienew = ienew + residueEnergy(savenrg,iprot2,ires2,nrgb,nrgf)
      DoneResNrg(iprot2,ires2) = .true.
      
      DoneResNrg(iprot1,ires1) = .false.
      DoneResNrg(iprot2,ires2) = .false.

c     call the montecarlo test...
      if (metropolis(ieold,ienew,movtyp)) then
c       if pass, then leave as it is
        titrationMove = .true.
        statstring = "PASS"
        if (updatenrg.eq.UPDATEENERGY) then
c         update the energies from residue ires2
          call updateResidueEnergies(iprot2)
c         recover the energy changes for residue ires1
          call popAllResidueEnergy
c         update the energies from residue ires1
          call updateResidueEnergies(iprot1)
c         update the extra energies
          call updateExtraEnergies
        endif
        call clearResidueEnergyStack
      else
c       if fail, then retun the residues to thier previous state...
        titrationMove = .false.
        statstring = "FAIL"
        if (convertRes(iprot1,ires1,ResNam(iprot1,ires1),oldres1)) then
          if (.not.convertRes(iprot2,ires2,ResNam(iprot2,ires2),oldres2)) then
            write(printstring,11) ResNam(iprot2,ires2),ResNum(iprot2,ires2)
            call printLine(WARNING,printstring)
          endif
        else
          write(printstring,10) ResNam(iprot1,ires1),ResNum(iprot1,ires1)
          call printLine(WARNING,printstring)
        endif
        call addAcid(iprot1,ires1)
        call takeBase(iprot1,ires1)
        call addBase(iprot2,ires2)
        call takeAcid(iprot2,ires2)
c       clear the residue energy buffer
        call clearResidueEnergyStack
      endif

      call clearMoveList()
      
      if (testenergy) call testEnergies      
      if (printMoveEvery.eq.0) return
                  
c     print the status of this move
      if (mod(movnum,printMoveEvery).eq.0) then
        write(printstring,12) movnum,iprot1,oldres1,ResNum(iprot1,ires1),iprot2,oldres2,
     .           ResNum(iprot2,ires2),ResNam(iprot1,ires1),ResNam(iprot2,ires2),
     .           OldTotalEnergy,OldTotalEnergyF,OldTotalEnergyB,statstring
12      format(i8,": TITRATION. Attempt",i2," ",a,"(",i3,") -> ",i2," ",
     .              a,"(",i3,") (now = ",a,"-",a,"), ",f12.3," | ",
     .              f12.3,f12.3," ",a4)
        call printLine(MOVE,printstring)
      endif
        
      return
      
      end
