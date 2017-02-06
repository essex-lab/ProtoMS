      
      logical function residueMove(updatenrg,movnum,movtyp)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'moves.inc'
      include 'energies.inc'
      include 'gb.inc'
      include 'quickpot.inc'
      include 'metropolis.inc'
c##########################################################
c
c     This function performs everything necessary for a
c     protein residue move... returns true if move is accepted
c
c     (C) Christopher Woods      5 February 2002
c
c##########################################################

      integer movnum,movtyp,updatenrg,savenrg
      double precision tempval
      integer iprot,ires,ityp
      double precision ieold,ienew,nrgb,nrgf
      double precision residueEnergy, ranRange
      double precision extraEnergy,inrgb,inrgf
      logical metropolis
      character*4 statstring
      integer ranInt,i,j
      double precision itmp

c     clear the move list
      call clearMoveList()      
            
c     randomly choose a protein..
      residueMove = .false.
      tempval = ZERO
      if (NProteins.le.0) return
 99   continue
      if (NProteins.eq.1) then
        iprot = 1
      else
        iprot = ranInt(1,NProteins)
      endif
      
c     randomly choose a residue...
      if (NResidues(iprot).le.0) then
        write(printstring,*) "No residues in protein ",iprot,"!"
        call printLine(WARNING,printstring)
        return
      endif

9     continue
c     This should make sure we don't pick fixed residues
      ires = ranInt(1,NResMove(iprot))
      ires = ResMove(iprot,ires)
c     make sure we don't try to move a fixed residue
      if (fixedRes(iprot,ires)) goto 9

      if (ProPrefSample) then
c        compute the owicki preferential sampling probability
c        for this residue, and compare with ran(0,1)
         itmp = ProProb(iprot,ires) / ProMaxProb
c        if we fail the test, then randomly pick another residue molecule
         if (itmp .le. ranRange(ZERO,ONE)) goto 99
      endif
      
c     increment the number of attempted moves for this residue
      NResTry(iprot,ires) = NResTry(iprot,ires) + 1

c     use a function to generate an integer between 1 and 3
c     1 = backbone move, 2 = sidechain move, 3 = both!
      if (FixedBac(iprot,ires)) then
        ityp = 2
      else
        ityp = ranInt(1,3)
      endif
      
c     add this residue to the move list
      call setResidueMoved(iprot,ires)      
c     may also need to add the next residue if this is not the last residue and
c     we will perform a backbone move
      if (ires.lt.NResidues(iprot).and.ityp.ne.2) call setResidueMoved(iprot,ires+1)
            
c     see whether we are updating the residue energies
      if (updatenrg.eq.UPDATEENERGY) then
        savenrg = OLD
      else
        savenrg = DONTCARE
      endif
      
c     get the old energy of this residue... - the internal angles of the
c     residue were calculated when it was built
      ieold = residueEnergy(savenrg,iprot,ires,nrgb,nrgf) + 
     .          ExtraEnergy(savenrg,inrgb,inrgf)
      DoneResNrg(iprot,ires) = .true.

c     if moving the backbone, then both this *and* the next residue 
c     must be rebuilt and evaluated.
      if (ires.lt.NResidues(iprot).and.ityp.ne.2) then
        call pushBackOldResidueEnergy
        ieold = ieold + residueEnergy(savenrg,iprot,ires+1,nrgb,nrgf)
      endif
      DoneResNrg(iprot,ires) = .false.

c     copy the old coordinates to storage...
      if (ires.lt.NResidues(iprot).and.ityp.ne.2)
     .                 call pushBackResidueCoords(iprot,ires+1)
      call pushBackResidueCoords(iprot,ires)
c     also need to create a residue - cut stack to save the residues that
c     were cut, in case of moving two residues at once...
      
      if (ityp.eq.1) then
c       move the backbone of the residue
        call moveAcid(iprot,ires,0,tempval)
      else if (ityp.eq.2) then
c       now move the sidechain
        call moveSideChain(iprot,ires)
      else
c       move both
        call moveAcid(iprot,ires,0,tempval)
        call moveSideChain(iprot,ires)
      endif
      
c     finally rebuild the rest of the residue from the zmatrix...
c     this recalculates the internal angles of the residue
      call buildResidue(iprot,ires)
      if (ires.lt.NResidues(iprot).and.ityp.ne.2)
     .                        call buildResidue(iprot,ires+1)

c     switch from old to new energy if necessary
      if (updatenrg.eq.UPDATEENERGY) savenrg = NEW
      
c     get the new energy
      ienew = residueEnergy(savenrg,iprot,ires,nrgb,nrgf) + 
     .             extraEnergy(savenrg,inrgb,inrgf)
      DoneResNrg(iprot,ires) = .true.
      if (ires.lt.NResidues(iprot).and.ityp.ne.2) then
        call pushBackNewResidueEnergy
        ienew = ienew + residueEnergy(savenrg,iprot,ires+1,nrgb,nrgf)
      endif
      DoneResNrg(iprot,ires) = .false.
      

c     if we are doing preferential sampling, then we must
c     use a biased MC test - this is based on the biased
c     probability of the new position, over the biased
c     probability of the old position
      if (ProPrefSample) then
        mt_bias = ProNewProb(iprot,ires) / ProProb(iprot,ires)
      endif
      
      
      if (metropolis(ieold,ienew,movtyp)) then
        residueMove = .true.
        statstring = "PASS"
        NResAcp(iprot,ires) = NResAcp(iprot,ires) + 1
c       update the energies...
        if (updatenrg.eq.UPDATEENERGY) then
          call updateResidueEnergies(iprot)
          if (ires.lt.NResidues(iprot).and.ityp.ne.2) then
            call popAllResidueEnergy
            call updateResidueEnergies(iprot)
          endif
          call updateExtraEnergies
        endif
c       we can leave the internal angles of the residue alone
c       we can now clear the residue coordinate stack
        call clearResidueCoordStack
        call clearResidueEnergyStack
        if (ProPrefSample) then
c         we need to save the new solvent probabilities!
           do i=1,NProteins
              do j=1,NResidues(i)
                 ProProb(i,j) = ProNewProb(i,j)
              enddo
          enddo
          ProMaxProb = ProNewMaxProb
          ProSumProb = ProNewSumProb
        endif
      else
        residueMove = .false.
        statstring = "FAIL"
c       copy back the old coordinates from the LIFO stack
c         the residue
        call popResidueCoords
c       JM this is annoying...
        if (doGB) call popPDACont
c         the next residue
        if (ires.lt.NResidues(iprot).and.ityp.ne.2) call popResidueCoords
c       we now need to update the residue internals to return them to
c       the pre-move state
        call updateResidueInternals(iprot,ires)
        if (ires.lt.NResidues(iprot).and.ityp.ne.2)
     .                call updateResidueInternals(iprot,ires+1)
        call clearResidueCoordStack
        call clearResidueEnergyStack
      endif

c     clear the move list if not doing a quick move
      if (.not. QuickMove) call clearMoveList()      
            
c     see whether we should test the energies...      
      if (testenergy) call testEnergies

      if (printMoveEvery.eq.0) return      
      if (mod(movnum,printMoveEvery).eq.0) then
        if (ityp.eq.1) then
          write(printstring,12) movnum,"BACKBONEMOVE ",iprot,ResNam(iprot,ires),
     .                                         ResNum(iprot,ires),oldTotalEnergy,statstring
        else if (ityp.eq.2) then
          write(printstring,12) movnum,"SIDECHAINMOVE",iprot,ResNam(iprot,ires),
     .                                         ResNum(iprot,ires),oldTotalEnergy,statstring
        else
          write(printstring,12) movnum,"RESIDUEMOVE  ",iprot,ResNam(iprot,ires),
     .                                         ResNum(iprot,ires),oldTotalEnergy,statstring
        endif        
12      format(i10,': ',a,1x,i3,1x,a,1x,i5,1x,g20.8,1x,a4)
        call printLine(MOVE,printstring)
        if (streamhandle(ENERGY).ge.0) call printEnergies(movnum)
      endif
      
      return
      end      
