      logical function volumeMove(updatenrg,movnum,movtyp)
      include 'dimensions.inc'
      include 'metropolis.inc'
      include 'printer.inc'
      include 'solvents.inc'
      include 'flags.inc'
      include 'moves.inc'
      include 'enums.inc'
      include 'energies.inc'
      include 'pressure.inc'
c###############################################################
c
c     This function performs everything necessary for a
c     volume move - thus allowing NPT simulations
c
c     (C) Christopher Woods    29 April 2003
c
c###############################################################

      integer updatenrg,movnum,savenrg
      integer movtyp,newtyp,i
      double precision oldvol,newvol,volchg
      double precision ienew,ieold,inrgb,inrgf
      double precision ranRange,totalEnergy,extraEnergy
      character*4 statstring
      logical metropolis

c     set everything as having moved
      call setAllMoved()      
      
      NVolTry = NVolTry + 1

c     get the current volume
      oldvol = volume
      
c     choose the new volume
      volchg = ranRange(-maxVolChange,maxVolChange)
      
      newvol = oldvol + volchg

c     see whether we are updating the solute energies
      if (updatenrg.eq.UPDATEENERGY) then
        savenrg = OLD
      else 
        savenrg = DONTCARE
      endif
            
c     everything will be moved
      call setAllMoved      
      
c     the energy before the move is the current total energy
c     stored in 'old'
      ieold = oldTotalEnergy

c     because we are not calculating the 'old' energy, set the old extra energies
c     to zero so that they are updated correctly - note that we are initialising
c     only the energies
      call initialiseExtraEnergies(.true.)
            
c     move the solvent molecules to generate the new volume...
      call moveVolume(newvol)

c     make the calculation go from old to new if necessary
      if (updatenrg.eq.UPDATEENERGY) savenrg = NEW
      
c     now recalculate the new energy...
      ienew = totalEnergy(savenrg,inrgb,inrgf) + extraEnergy(savenrg,inrgb,inrgf)

c     make sure that we use NPT acceptance test
      if (movtyp.le.0) then
        newtyp = NPTMCTEST
      else
        newtyp = movtyp
      endif

      mt_vold = oldvol
      mt_vnew = newvol

      if (metropolis(ieold,ienew,newtyp)) then
        volumemove = .true.
        statstring = "PASS"
        NVolAcp = NVolAcp + 1
        if (updatenrg.eq.UPDATEENERGY) then
          call updateTotalEnergies
          call updateExtraEnergies
        endif
        if (PrefSample) then
c         we need to save the new solvent probabilities!
          do i=1,NSolvents
            SvnProb(i) = SvnNewProb(i)
          enddo
          SvnMaxProb = SvnNewMaxProb
          SvnSumProb = SvnNewSumProb
        endif
      else
c       change the volume back again!
        call reverseVolumeMove
        volumemove = .false.
        statstring = "FAIL"
      endif

c     clear the move list
      call clearMoveList
      
      if (testenergy) call testEnergies

c     print the status of this move
      if (printMoveEvery.eq.0) return
      if (mod(movnum,printMoveEvery).eq.0) then                  
        write(printstring,12) movnum,newvol,oldTotalEnergy,statstring
        call printLine(MOVE,printstring)
12      format(i10,": VOLUMEMOVE ",f15.5," A^3 ",g20.8,1x,a4)

        if (mod(movnum,printMoveEvery).eq.0) call printEnergies(movnum)
      endif

      return
      end
