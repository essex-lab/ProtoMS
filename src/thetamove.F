      logical function thetaMove(updatenrg,movnum,movtyp)
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
      include 'simparams.inc'
c###############################################################
c
c     This function performs everything necessary for a 
c     successful GC-solute theta move. 
c
c     (C) Michael Bodnarchuk   22nd October 2009
c###############################################################

      integer updatenrg,movnum,movtyp,igcs,ranInt,savenrg
      double precision GCSoluteEnergy,cnrgf,cnrgb,ljnrgf,ljnrgb,ranRange
      double precision ieold,ienew,thetaold,nrg,ranval
      character*4 statstring
      logical metropolis
 
      thetaMove = .false.

c     Select a test particle at random
      igcs = ranInt(1,NGCSolutes)
      thetaold = thetaparam(igcs)
 
      if (updatenrg.eq.UPDATEENERGY) then
        savenrg = OLD
      else
        savenrg = DONTCARE
      endif
      
      NTheTry(igcs) = NTheTry(igcs) + 1
c     Get the old energy of the theta molecule 
      ieold = GCsoluteEnergy(savenrg,igcs,cnrgf,cnrgb,ljnrgf,ljnrgb)
c     Change the theta value of the molecule according to the random number
10    ranval = ranRange(-GCSADel(igcs),GCSADel(igcs))
c     Note that this uses ranRange, so it adds or removes +/- 0.15
      thetaparam(igcs) = thetaold + ranval
      if (thetaparam(igcs).gt.1) then
      goto 10
      endif
      if (thetaparam(igcs).lt.0) then 
      goto 10
      endif
      
c     See if the new configurational energy needs to be updated...
      if (updatenrg.eq.UPDATEENERGY) savenrg = NEW      
      ienew = GCsoluteEnergy(savenrg,igcs,cnrgf,cnrgb,ljnrgf,ljnrgb)

      if (metropolis(ieold,ienew,movtyp)) then
        thetaMove = .true.
        statstring = "PASS"          
        nrg = ienew  
        NTheAcp(igcs) = NTheAcp(igcs) + 1
c       update the energies
        if (updatenrg.eq.UPDATEENERGY) then
           call updateGCSoluteEnergies(igcs)
        endif
      else
        thetaMove = .false.
        statstring = "FAIL"
        nrg = ieold
        thetaparam(igcs) = thetaold
      endif

c     savepoint dictates whether we save or not...
      if (thetaparam(igcs).gt.thres) then
         call savepoint(igcs,movnum)
      endif

      if (printMoveEvery.eq.0) return
      if (mod(movnum,printMoveEvery).eq.0) then
        write(printstring,12) movnum,igcs,GCSName(igcs),oldTotalEnergy,statstring,thetaparam(igcs)
12      format(i10,": THETAMOVE ",i10," ",a," ",g20.8,1x,a4,g20.8)
        call printLine(MOVE,printstring)
        if (streamhandle(ENERGY).ge.0) call printEnergies(movnum)
      endif
      
      return
      end
