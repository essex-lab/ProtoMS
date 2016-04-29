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

      integer updatenrg,movnum,movtyp,igcs,ranInt,savenrg
      double precision GCSoluteEnergy,cnrgf,cnrgb,ljnrgf,ljnrgb
      double precision ieold,ienew,thetaold,nrg
      character*4 statstring
      logical metropolis2

      
      insertionMove = .false.

c     Select a test particle at random
10    igcs = ranInt(1,NGCSolutes)
c     Check to see that it's not already on
      if (thetaparam(igcs).eq.1) then
c      PRINT*, "oops",movnum
      return
      endif
      thetaold = thetaparam(igcs)
 
      if (updatenrg.eq.UPDATEENERGY) then
        savenrg = OLD
      else
        savenrg = DONTCARE
      endif
      
      NInsTry(igcs) = NInsTry(igcs) + 1
c     Get the old energy of the theta molecule - should be zero 
      ieold = GCsoluteEnergy(savenrg,igcs,cnrgf,cnrgb,ljnrgf,ljnrgb)
      if (ieold.gt.0) then
      print*, "This isn't good!"
      endif
c     Change the theta value of the molecule to zero
      thetaparam(igcs) = 1
      call pushBackGCSoluteCoords(igcs)



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
        call clearGCSoluteCoordStack
      else
        insertionMove = .false.
        statstring = "FAIL"
        nrg = ieold
        thetaparam(igcs) = thetaold
        GCSCavity = GCSCavity - 1
        call popGCSoluteCoords
        if (doGB) call popPDACont
        call clearGCSoluteCoordStack
        
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
