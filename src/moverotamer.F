
      subroutine moveRotamer(iprot,ires)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'rotamer.inc'
      include 'printer.inc'
c####################################################
c
c     This subroutine actually performs the 
c     rotamer library move on residue ires of 
c     protein iprot.
c
c     (C) Christopher Woods, October 2005
c
c####################################################

      integer iprot,ires
      integer nsc
      integer j,irot
      integer itmpl
      double precision ranval, ranRange, grnd
      integer ranInt
      logical dorotamer

c     check that this residue exists
      if (iprot.le.0 .or. iprot.gt.NProteins) return
      if (ires.le.0 .or. ires.gt.NResidues(iprot)) return

c     there are two types of rotamer move. The first is where the dihedrals
c     are changed to all lie within 10 degrees of a rotamer state.

c     The second is where all of the dihedrals are changed by large amounts
c     (by +/- 180 degrees)

      itmpl = ResRotamerID(iprot,ires)

c     first, see if rotamer library moves are possible for this residue
      if (itmpl .ne. 0 .and. BigMoveProb .lt. ONE) then
c       there are some rotamer dihedrals and states for this residue!
        
c       see if we should do a rotamer move, or do a big dihedral move
        dorotamer = grnd() .ge. BigMoveProb
      else
        if (BigMoveProb .eq. ZERO) then
          write(printstring,*) "No rotamer library for residue ",ResNam(iprot,ires),
     .         " and BigMoves are disabled. Skipping the rotamer move!"
          call printLine(WARNING,printstring)
          return
        endif
      
        dorotamer = .false.
      endif

      if (dorotamer) then
c       choose a rotamer state to target
        irot = ranInt(1, NRotamerStates(itmpl))

c       move the dihedrals in the side chain to be near the rotamer values of that state
        nsc = NSCAtoms(iprot,ires)
        if (nsc.gt.0) then
          do j=1,nsc
            if (ZMovSC(iprot,ires,j,3) .and. RotamerID(iprot,ires,j) .ne. 0) then
c             this is a rotamer dihedral - set it to within +/- RotamerDelta of
c             the rotamer value
              ranval = ranRange(-RotamerDelta, RotamerDelta)
              ZValSC(iprot,ires,j,3) = RotamerVal(itmpl,irot,RotamerID(iprot,ires,j))
     .                                   + ranval
            endif
          enddo
        endif

c       done!
      else
c       do the big dihedral moves
        nsc = NSCAtoms(iprot,ires)
        if (nsc.gt.0) then
          do j=1,nsc
            if (ZMovSC(iprot,ires,j,3)) then
              ranval = ranRange(-BigMoveDelta, BigMoveDelta)
              ZValSC(iprot,ires,j,3) = ZValSC(iprot,ires,j,3) + ranval
            endif
          enddo
        endif
      endif     

      return
      end
