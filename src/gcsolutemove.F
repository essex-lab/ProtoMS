      
      logical function gcsoluteMove(updatenrg,movnum,movtyp)
      include 'dimensions.inc'
      include 'metropolis.inc'
      include 'printer.inc'
      include 'constants.inc'
      include 'solvents.inc'
      include 'solutes.inc'
      include 'flags.inc'
      include 'moves.inc'
      include 'enums.inc'
      include 'energies.inc'
      include 'gb.inc'
      include 'quickpot.inc'
      include 'stacks.inc'
c##################################################################
c
c     This function performs everything necessary for a GC-solute
c     move.. Returns true if the move was accepted
c     This is a SPECIALISED GCMC FUNCTION
c     MSB Mar 2009
cdsasa
c##################################################################

      integer updatenrg,movnum,movtyp,savenrg
      integer ind,count,i,igcs
      double precision ieold,blah,cnrgf,cnrgb,ljnrgf,ljnrgb,autofail
      double precision GCSoluteEnergy
      double precision ienew
      double precision nrg
      logical metropolis
      character*4 statstring
      integer ranInt

c     clear the move list
      call clearMoveList()      
      
c     randomly choose a solvent
      gcsoluteMove = .false.

c     different selection selection logic if we are doing jaws or gcmc
      if (doJaws) then
         igcs = ranInt(1,NGCSolutes)
      else
         if (GCSCavity.eq.0) return

c        Select from inserted gcsolutes randomly
         ind = ranInt(1,GCSCavity)

c        Now find the ind'th on gcsolute
         igcs = 0
         count = 0
         do i=1,NGCSolutes
            if (thetaparam(i).eq.1) then
               count = count + 1
            endif
            if (count.eq.ind) then
               igcs = i
               goto 10
            endif
         enddo
      endif

      
c     add this solute to the move list
 10   call setGCSoluteMoved(igcs)      
c     see whether we are updating the energies
      if (updatenrg.eq.UPDATEENERGY) then
        savenrg = OLD
      else
        savenrg = DONTCARE
      endif
     
c      do p=1,NGCSolutes
c      lah = GCsoluteEnergy(savenrg,p,cnrgf,cnrgb,ljnrgf,ljnrgb)
c      PRINT*, "p",lah
c      enddo
c     get the old energy of the solute
c      ieold = totalEnergy(savenrg,nrgb,nrgf)  
        ieold = GCsoluteEnergy(savenrg,igcs,cnrgf,cnrgb,ljnrgf,ljnrgb)
c     copy the coordinates to storage
      call pushBackGCSoluteCoords(igcs)

c     increment the number of attempted moves of the solvent
      NGCSTry(igcs) = NGCSTry(igcs) + 1      
      
c     move the solvent
      call moveGCSolute(igcs,ieold,ienew)
c     change to NEW energy if necessary
      if (updatenrg.eq.UPDATEENERGY) savenrg = NEW
      
c     get the new energy
        blah = GCsoluteEnergy(savenrg,igcs,cnrgf,cnrgb,ljnrgf,ljnrgb)
c        PRINT*, "blah2" ,blah
        blah = ienew + blah + autofail
     
      if (metropolis(ieold,blah,movtyp)) then
        gcsoluteMove = .true.
        statstring = "PASS"      
        nrg = blah      
        NGCSAcp(igcs) = NGCSAcp(igcs) + 1
c       update the energies
        if (updatenrg.eq.UPDATEENERGY) then
          call updateGCSoluteEnergies(igcs)
        endif

c       now clear the solvent coordinate stack
        call clearGCSoluteCoordStack
      else
        gcsoluteMove = .false.
        statstring = "FAIL"  
        nrg = ieold    
c       copy back the old coords from the stack
        call popGCSoluteCoords
c       JM this is annoying...
        if (doGB) call popPDACont
        call clearGCSoluteCoordStack
      endif

c     clear the move lists
      if (.not. QuickMove) call clearMoveList()      
      
      if (testenergy) call testEnergies
c         SvnCavity = 0
c         SvnReservoir = 0
c     print the status of this move
      if (printMoveEvery.eq.0) return
      if (mod(movnum,printMoveEvery).eq.0) then
        write(printstring,12) movnum,igcs,GCSName(igcs),nrg,statstring
12      format(i10,": GCSOLUTEMOVE ",i10," ",a," ",g20.8,1x,a4)
        call printLine(MOVE,printstring)
        if (streamhandle(ENERGY).ge.0) call printEnergies(movnum)
      endif
      
      return
      end
