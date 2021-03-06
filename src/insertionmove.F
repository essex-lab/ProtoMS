      logical function insertionMove(updatenrg,movnum,movtyp)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'stacks.inc'
      include 'solvents.inc'
      include 'metropolis.inc'
      include 'constants.inc'
      include 'solutes.inc'
      include 'flags.inc'
      include 'moves.inc'
      include 'enums.inc'
      include 'energies.inc'
      include 'gb.inc'
      include 'quickpot.inc'


c###############################################################
c
c     This function performs everything necessary for a 
c     successful GC insertion move. Similar to the JAWS routine, 
c     it switches the theta value of a molecule from 0 to 1
c
c     (C) Michael Bodnarchuk   15th September 2011
c
c###############################################################

      integer updatenrg,movnum,movtyp,ind,i,count,igcs,ranInt,savenrg
      double precision GCSoluteEnergy,cnrgf,cnrgb,ljnrgf,ljnrgb
      double precision zeroGCSoluteEnergy
      double precision ieold,ienew,thetaold,nrg
      character*4 statstring
      logical metropolis2 
      
      insertionMove = .false.

      if (GCSCavity.eq.NGCSolutes) then
         write(printstring,*) 'No spare gcsolutes to insert. This will bias your simulation results.'
         call PrintLine(WARNING,printstring)
         return
      endif

c     for purposes of insertion all off GCS are interchangeable
c     so just find the first one that is switched off
      do i=1,NGCSolutes
         if (thetaparam(i).eq.0) then
            igcs = i
            goto 10
         endif
      enddo

 10   thetaold = thetaparam(igcs)
 
      if (updatenrg.eq.UPDATEENERGY) then
        savenrg = OLD
      else
        savenrg = DONTCARE
      endif
c      call scramblewater(igcs)
      
      NInsTry(igcs) = NInsTry(igcs) + 1
c     Get the old energy of the theta molecule - should be zero 
c      ieold = GCsoluteEnergy(savenrg,igcs,cnrgf,cnrgb,ljnrgf,ljnrgb)
c$$$      ieold = ZERO
c$$$      cnrgf = ZERO
c$$$      ljnrgf = ZERO
c$$$      ljnrgb = ZERO
      ieold = zeroGCSoluteEnergy(savenrg,igcs,cnrgf,cnrgb,ljnrgf,ljnrgb)
      if (ieold.gt.0) then
      print*, "This isn't good!"
      endif
c     Change the theta value of the molecule to zero
      thetaparam(igcs) = 1
        
      call scramblewater(igcs)

      GCSCavity = GCSCavity + 1
c     See if the new configurational energy needs to be updated...
      if (updatenrg.eq.UPDATEENERGY) savenrg = NEW      
      ienew = GCsoluteEnergy(savenrg,igcs,cnrgf,cnrgb,ljnrgf,ljnrgb)

      if (metropolis2(ieold,ienew,movtyp)) then
        insertionMove = .true.
        statstring = "PASS"          
        nrg = ienew  
        NInsAcp(igcs) = NInsAcp(igcs) + 1
        
c       update the energies
        if (updatenrg.eq.UPDATEENERGY) then
           call updateGCSoluteEnergies(igcs)
        endif
c        if (thetaparam(igcs).gt.0.95) then
c           call savepoint(igcs,movnum)
c        endif
      else
        insertionMove = .false.
        statstring = "FAIL"
        nrg = ieold
        thetaparam(igcs) = thetaold
        GCSCavity = GCSCavity - 1
        if (doGB) call popPDACont
        
      endif
c       if (thetaparam(igcs).gt.0.950) then
c           call savepoint(igcs,movnum)
c        endif
       if (printMoveEvery.eq.0) return
      if (mod(movnum,printMoveEvery).eq.0) then
        write(printstring,12) movnum,igcs,GCSName(igcs),oldTotalEnergy,statstring,thetaparam(igcs)
12      format(i10,": INSERTIONMOVE ",i10," ",a," ",g20.8,1x,a4,g20.8)
        call printLine(MOVE,printstring)
        if (streamhandle(ENERGY).ge.0) call printEnergies(movnum)
      endif
      return
      end
