      logical function sampleMove(updatenrg,movnum,movtyp)
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
c     successful theta-sampling move. 
c     A random change is made, and a bias potential applied. If the
c     accepted value is greater than 0.95 or less than 0.05 then this is
c     held for analysis.
c
c     (C) Michael Bodnarchuk   5th January 2010
c
c###############################################################

      integer updatenrg,movnum,movtyp,igcs,ranInt,savenrg
      double precision GCSoluteEnergy,cnrgf,cnrgb,ljnrgf,ljnrgb,ranval,ranRange
      double precision ieold,ienew,thetaold,nrg,potential,biasPotential
      character*4 statstring
      logical metropolis
      
      sampleMove = .false.

c     Select a test particle at random
      igcs = ranInt(1,NGCSolutes)
      thetaold = thetaparam(igcs)

      if (updatenrg.eq.UPDATEENERGY) then
        savenrg = OLD
      else
        savenrg = DONTCARE
      endif
      NSamTry(igcs) = NSamTry(igcs) + 1

c     Get the old energy of the theta molecule 
      ieold = GCsoluteEnergy(savenrg,igcs,cnrgf,cnrgb,ljnrgf,ljnrgb)

c     Get the bias
      potential = biasPotential(igcs,potential)

c     Add the bias onto the old energy of the gcsolute
      ieold = ieold + potential

c     Change the theta value of the molecule according to the random number

10    ranval = ranRange(-GCSADel(igcs),GCSADel(igcs))
      thetaparam(igcs) = thetaold + ranval
      if (thetaparam(igcs).gt.1) then
      goto 10
      endif
      if (thetaparam(igcs).lt.0) then
      goto 10
      endif

c     See if the new configurational energy needs to be updated...
      if (updatenrg.eq.UPDATEENERGY) savenrg = NEW      

c     Get the new energy
      ienew = GCsoluteEnergy(savenrg,igcs,cnrgf,cnrgb,ljnrgf,ljnrgb)

c     Get the bias
      potential = biasPotential(igcs,potential)

c     Add the bias onto the new energy of the gcsolute
      ienew = ienew + potential

      if (metropolis(ieold,ienew,movtyp)) then
        sampleMove = .true.
        statstring = "PASS"          
        nrg = ienew  
        NSamAcp(igcs) = NSamAcp(igcs) + 1
c       update the energies
        if (updatenrg.eq.UPDATEENERGY) then
           call updateGCSoluteEnergies(igcs)
        endif

c       Track to see if the theta value is high/low 
        if (thetaparam(igcs).ge.0.95) then
           call savehigh(igcs,movnum)
        endif
        if (thetaparam(igcs).le.0.05) then
           call savelow(igcs,movnum)
        endif

      else
        sampleMove = .false.
        statstring = "FAIL"
        nrg = ieold
        thetaparam(igcs) = thetaold

c       Track to see if the theta value is high/low
        if (thetaparam(igcs).ge.0.95) then
           call savehigh(igcs,movnum)
        endif
        if (thetaparam(igcs).le.0.05) then
           call savelow(igcs,movnum)
        endif
      endif

      if (printMoveEvery.eq.0) return
      if (mod(movnum,printMoveEvery).eq.0) then
        write(printstring,12) movnum,igcs,GCSName(igcs),oldTotalEnergy,statstring,thetaparam(igcs)
12      format(i10,": SAMPLEMOVE ",i10," ",a," ",g20.8,1x,a4,g20.8)
        call printLine(MOVE,printstring)
        if (streamhandle(ENERGY).ge.0) call printEnergies(movnum)
      endif
      
      return
      end
