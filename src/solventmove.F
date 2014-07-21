      
      logical function solventMove(updatenrg,movnum,movtyp)
      include 'dimensions.inc'
      include 'metropolis.inc'
      include 'printer.inc'
      include 'constants.inc'
      include 'solvents.inc'
      include 'flags.inc'
      include 'moves.inc'
      include 'enums.inc'
      include 'energies.inc'
      include 'gb.inc'
      include 'quickpot.inc'
      include 'stacks.inc'
c##################################################################
c
c     This function performs everything necessary for a solvent
c     move.. Returns true if the move was accepted
c     This is a SPECIALISED GCMC FUNCTION
c     MSB Mar 2009
cdsasa
c##################################################################

      integer updatenrg,movnum,movtyp,savenrg
      integer isvn
      double precision ieold,itmp
      double precision solventEnergy,extraEnergy
      double precision ranRange,ienew
      double precision nrgf,nrgb,nrg
      logical metropolis
      character*4 statstring
      integer ranInt,i

c     clear the move list
      call clearMoveList()      
      
c     randomly choose a solvent
      solventMove = .false.
      
      if (NSolvents.eq.1) then
        isvn = 1
      else
99      continue
        isvn = ranInt(1,NSolvents)

        if (PrefSample) then
c         compute the owicki preferential sampling probability
c         for this solvent, and compare with ran(0,1)
          itmp = SvnProb(isvn) / SvnMaxProb
c         if we fail the test, then randomly pick another solvent molecule
          if (itmp .le. ranRange(ZERO,ONE)) goto 99
        endif
      endif

c     add this solvent to the move list
      call setSolventMoved(isvn)      
            
c     see whether we are updating the energies
      if (updatenrg.eq.UPDATEENERGY) then
        savenrg = OLD
      else
        savenrg = DONTCARE
      endif
      
c     get the old energy of the solvent
      ieold = solventEnergy(savenrg,isvn,nrgb,nrgf)  + 
     .          extraEnergy(savenrg,nrgb,nrgf)
               
c     copy the coordinates to storage
      call pushBackSolventCoords(isvn)

c     increment the number of attempted moves of the solvent
      NSvnTry(isvn) = NSvnTry(isvn) + 1      
      
c     move the solvent
      call moveSolvent(isvn,ienew)
      
c     change to NEW energy if necessary
      if (updatenrg.eq.UPDATEENERGY) savenrg = NEW
      
c     get the new energy
      ienew = ienew + solventEnergy(savenrg,isvn,nrgb,nrgf)  + 
     .          extraEnergy(savenrg,nrgb,nrgf)
                

c     if we are doing preferential sampling, then we must
c     use a biased MC test - this is based on the biased
c     probability of the new position, over the biased
c     probability of the old position
      if (PrefSample) then
        mt_bias = SvnNewProb(isvn) / SvnProb(isvn)
      endif
      
c     now perform the metropolis test
c      write(printstring,*) 'ieold ',ieold,' ienew ',ienew
c      call printLine(DEBUG,printstring)
      if (metropolis(ieold,ienew,movtyp)) then
        solventMove = .true.
        statstring = "PASS"      
        nrg = ienew      
        NSvnAcp(isvn) = NSvnAcp(isvn) + 1
c       update the energies
        if (updatenrg.eq.UPDATEENERGY) then
          call updateSolventEnergies  
          call updateExtraEnergies

        endif

c       now clear the solvent coordinate stack
        call clearSolventCoordStack
        if (PrefSample) then
c         we need to save the new solvent probabilities!
          do i=1,NSolvents
            SvnProb(i) = SvnNewProb(i)
          enddo
          SvnMaxProb = SvnNewMaxProb
          SvnSumProb = SvnNewSumProb
        endif
      else
        solventMove = .false.
        statstring = "FAIL"  
        nrg = ieold    
c       copy back the old coords from the stack
        call popSolventCoords
c       JM this is annoying...
        if (doGB) call popPDACont
        call clearSolventCoordStack
      endif

c     clear the move lists
      if (.not. QuickMove) call clearMoveList()      
      
      if (testenergy) call testEnergies

c     print the status of this move
      if (printMoveEvery.eq.0) return
      if (mod(movnum,printMoveEvery).eq.0) then
        write(printstring,12) movnum,isvn,SvnName(isvn),nrg,statstring
12      format(i10,": SOLVENTMOVE ",i10," ",a," ",g20.8,1x,a4)
        call printLine(MOVE,printstring)
        if (streamhandle(ENERGY).ge.0) call printEnergies(movnum)
      endif
      
      return
      end
