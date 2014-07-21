      double precision function multipleMovesEnergy(savenrg,nrgb,nrgf)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'solvents.inc'
      include 'movelist.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'energies.inc'
c#######################################
c  This routine returns the sum of 
c  the energies of solutes,residues,
c  solvents...that have moved
c
c
c     (c) Julien Michel April 2005
c######################################
      integer i,j,nres,savenrg
      double precision inrg,inrgb,inrgf
      double precision nrgb,nrgf
      double precision residueEnergy
      double precision soluteEnergy
      double precision solventEnergy
      double precision extraEnergy
      logical hasResidueMoved
      logical hasSoluteMoved
      logical hasSolventMoved

      multipleMovesEnergy = ZERO
      nrgb = ZERO
      nrgf = ZERO
      inrg = ZERO
      inrgf = ZERO
      inrgb = ZERO

c     Make sure all FEP solutes are synchronised...
      do i=1,NSolutes
        if (IgnoreSolute(i)) continue
        if (DifferentLambdas .and. NSolVariable(i).gt.0) then
          call synchFEPSolutes(i)
        endif
      enddo

      if (NProteins.gt.0) then
c       loop over all proteins and get them to return thier
c       intramolecular energy...
        do i=1,NProteins
c         loop over all residues in the protein and see if they have moved
          nres = NResidues(i)
          do j=1,nres
            if (hasResidueMoved(i,j)) then 
              inrg = residueEnergy(savenrg,i,j,inrgb,inrgf)
              DoneResNrg(i,j) = .true.
              multipleMovesEnergy = multipleMovesEnergy + inrg
              nrgb = nrgb + inrgb
              nrgf = nrgf + inrgf
c              write(printstring,*) 'SaveNrg ',savenrg,' residue ',i,'-',j,
c     .          ' inrg ',inrg, ' f ',nrgf,' b ',inrgb
c              call printLine(DEBUG,printstring)
              if (savenrg .eq. NEW) then 
                call pushBackNewResidueEnergy
              else if (savenrg .eq. OLD) then
                call pushBackOldResidueEnergy
              endif
            endif
          enddo
        enddo
      endif

c     now do the same for all of the solutes...
      if (NSolutes.gt.0) then
        do i=1,NSolutes
          if (.not.IgnoreSolute(i)) then
c           see if the solute has moved
            if (hasSoluteMoved(i)) then
              inrg = soluteEnergy(savenrg,i,inrgb,inrgf)
              DoneSolNrg(i) = .true.
c              write(printstring,*) 'SaveNrg ',savenrg,' solute ',i,
c     .          ' inrg ',inrg, ' f ',inrgf,' b ',inrgb
c              call printLine(DEBUG,printstring)
              multipleMovesEnergy = multipleMovesEnergy + inrg
              nrgb = nrgb + inrgb
              nrgf = nrgf + inrgf
            endif
          endif        
        enddo
      endif

c     now do the same for all the solvents...
      if (NSolvents.gt.0) then
        do i=1,NSolvents
c           see if the solute has moved
          if (hasSolventMoved(i)) then
            inrg = solventEnergy(savenrg,i,inrgb,inrgf)
            DoneSvnNrg(i) = .true.
            multipleMovesEnergy = multipleMovesEnergy + inrg
            nrgb = nrgb + inrgb
            nrgf = nrgf + inrgf
          endif
        enddo
      endif

c     Now do the extra energies
      inrg = extraEnergy(savenrg,inrgb,inrgf)
c      write(printstring,*) 'The ExtraEnergy ',inrg
c      call printLine(DEBUG,printstring)
      multipleMovesEnergy = multipleMovesEnergy + inrg
      nrgb = nrgb + inrgb
      nrgf = nrgf + inrgf      

c     Clear the DoneFlags
      if (NProteins.gt.0) then
c       loop over all proteins and get them to return thier
c       intramolecular energy...
        do i=1,NProteins
          nres = NResidues(i)
          do j=1,nres
            if (hasResidueMoved(i,j)) then 
              DoneResNrg(i,j) = .false.
            endif
          enddo
        enddo
      endif      
      if (NSolutes.gt.0) then
        do i=1,NSolutes
          if (.not.IgnoreSolute(i)) then
c           see if the solute has moved
            if (hasSoluteMoved(i)) then
              DoneSolNrg(i) = .false.
            endif
          endif        
        enddo
      endif      
c     now do the same for all the solvents...
      if (NSolvents.gt.0) then
        do i=1,NSolvents
c           see if the solute has moved
          if (hasSolventMoved(i)) then
            DoneSvnNrg(i) = .false.
          endif
        enddo
      endif

      return
      end
