
      logical function rotamerMove(updatenrg,movnum,movtyp)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'moves.inc'
      include 'energies.inc'
c###########################################################
c
c     This function is used to perform rotamer library
c     moves on protein residues. Note that these moves
c     do not respect detailed balance, so are not suitable
c     for equilibrium simulations (e.g. Monte Carlo free
c     energy calculations). These moves should only
c     be used if you are relaxing a system, or if 
c     you are doing docking.
c
c     (C) Christopher Woods,  October 2005
c
c###########################################################

      integer updatenrg,movnum,movtyp
      integer iprot,ires
      integer i,savenrg
      double precision ieold,ienew,nrgb,nrgf,inrgb,inrgf
      integer ranInt
      double precision residueEnergy,extraEnergy
      logical haveflex,metropolis
      character*4 statstring

c     clear the move list
      call clearMoveList()

c     randomly choose a protein..
      rotamerMove = .false.

      if (NProteins.le.0) return

c     make sure that there are some flexible residues, or this will go into 
c     an infinite loop!
      haveflex = .false.
      do i=1,NProteins
        if (NResMove(i).gt.0) haveflex = .true.
      enddo

      if (.not. haveflex) then
c       we have a problem - there are no flexible protein residues!               
        call printLine(WARNING,"Attempting a residue move when there are no flexible residues!")
        return
      endif

8     continue
      if (NProteins.eq.1) then
        iprot = 1
        if (NResMove(iprot) .le. 0) return
      else
        iprot = ranInt(1,NProteins)
      endif
      
c     skip this protein if there are no flexible residues
      if (NResMove(iprot) .le. 0) goto 8

c     randomly choose a residue...
      if (NResidues(iprot).le.0) goto 8

c     This should make sure we don't pick fixed residues
      ires = ranInt(1,NResMove(iprot))
      ires = ResMove(iprot,ires)

c     make sure we don't try to move a fixed residue
      if (fixedRes(iprot,ires)) then
        write(printstring,*) "Residue ",iprot,"-",ires," is fixed, but is listed ",
     .                       "as a flexible residue. This is a program bug! Skipping the move!"
        call printLine(WARNING,printstring)
        return
      endif

c     increment the number of attempted moves for this residue
      NResTry(iprot,ires) = NResTry(iprot,ires) + 1
      
c     add this residue to the move list
      call setResidueMoved(iprot,ires)      

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

c     copy the old coordinates to storage...
      call pushBackResidueCoords(iprot,ires)

c     perform the rotamer move
      call moveRotamer(iprot,ires)
      
c     rebuild the residue
      call buildResidue(iprot,ires)

c     switch from old to new energy if necessary
      if (updatenrg.eq.UPDATEENERGY) savenrg = NEW

c     get the new energy
      ienew = residueEnergy(savenrg,iprot,ires,nrgb,nrgf) + 
     .             extraEnergy(savenrg,inrgb,inrgf)

      if (metropolis(ieold,ienew,movtyp)) then
        rotamerMove = .true.
        statstring = "PASS"
        NResAcp(iprot,ires) = NResAcp(iprot,ires) + 1

c       update the energies...
        if (updatenrg.eq.UPDATEENERGY) then
          call updateResidueEnergies(iprot)
          call updateExtraEnergies
        endif

c       we can leave the internal angles of the residue alone
c       we can now clear the residue coordinate stack
        call clearResidueCoordStack
        call clearResidueEnergyStack
      else
        rotamerMove = .false.
        statstring = "FAIL"
c       copy back the old coordinates from the LIFO stack
        call popResidueCoords

c       we now need to update the residue internals to return them to
c       the pre-move state
        call updateResidueInternals(iprot,ires)

        call clearResidueCoordStack
        call clearResidueEnergyStack
      endif

c     clear the move list
      call clearMoveList()      
            
c     see whether we should test the energies...      
      if (testenergy) call testEnergies

      if (printMoveEvery.eq.0) return      
      if (mod(movnum,printMoveEvery).eq.0) then
        write(printstring,12) movnum,"ROTAMERMOVE  ",iprot,ResNam(iprot,ires),
     .                        ResNum(iprot,ires),oldTotalEnergy,statstring
12      format(i10,': ',a,1x,i3,1x,a,1x,i5,1x,g20.8,1x,a4)
        call printLine(MOVE,printstring)
        if (streamhandle(ENERGY).ge.0) call printEnergies(movnum)
      endif
      
      return
      end      
