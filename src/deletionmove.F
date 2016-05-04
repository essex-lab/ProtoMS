      logical function deletionMove(updatenrg,movnum,movtyp)
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
c     successful GC deletion move. Similar to the JAWS routine, 
c     it switches the theta value of a molecule from 1 to 0
c
c     (C) Michael Bodnarchuk   15th September 2011
c
c###############################################################

      integer updatenrg,movnum,movtyp,ind,count,i,igcs,ranInt,savenrg
      double precision GCSoluteEnergy,cnrgf,cnrgb,ljnrgf,ljnrgb
      double precision ieold,ienew,thetaold,nrg
      character*4 statstring
      logical metropolis3

      
      deletionMove = .false.
      if (GCSCavity.eq.0) return

c     Select from inserted gcsolutes randomly
      ind = ranInt(1,GCSCavity)

c     Now find the ind'th on gcsolute
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
      
 10   thetaold = thetaparam(igcs)
 
      if (updatenrg.eq.UPDATEENERGY) then
        savenrg = OLD
      else
        savenrg = DONTCARE
      endif
      
      
      NDelTry(igcs) = NDelTry(igcs) + 1
c     Get the old energy of the theta molecule - should be zero 
      ieold = GCsoluteEnergy(savenrg,igcs,cnrgf,cnrgb,ljnrgf,ljnrgb)

      if (ieold.eq.0) then
      write(printstring,*) "Energy of gcsolute isn't >zero! BUG"
      call printLine(WARNING,printstring)
      endif

c     Change the theta value of the molecule to zero
      thetaparam(igcs) = 0
c     See if the new configurational energy needs to be updated...
      if (updatenrg.eq.UPDATEENERGY) savenrg = NEW      
      ienew = GCsoluteEnergy(savenrg,igcs,cnrgf,cnrgb,ljnrgf,ljnrgb)

       if (metropolis3(ieold,ienew,movtyp)) then
        deletionMove = .true.
        statstring = "PASS"          
        nrg = ienew  
        NDelAcp(igcs) = NDelAcp(igcs) + 1
        GCSCavity = GCSCavity - 1
c       update the energies
        if (updatenrg.eq.UPDATEENERGY) then
           call updateGCSoluteEnergies(igcs)
        endif
      else
        deletionMove = .false.
        statstring = "FAIL"
        nrg = ieold
        thetaparam(igcs) = thetaold
      endif
c      PRINT*, "carried GCS",GCSCavity
      if (printMoveEvery.eq.0) return
      if (mod(movnum,printMoveEvery).eq.0) then
        write(printstring,12) movnum,igcs,GCSName(igcs),oldTotalEnergy,statstring,thetaparam(igcs)
12      format(i10,": DELETIONMOVE ",i10," ",a," ",g20.8,1x,a4,g20.8)
        call printLine(MOVE,printstring)
        if (streamhandle(ENERGY).ge.0) call printEnergies(movnum)
      endif
      
      return
      end
