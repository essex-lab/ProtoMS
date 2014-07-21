      
      logical function soluteMove(updatenrg,movnum,movtyp)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'solutes.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'moves.inc'
      include 'energies.inc'
      include 'gb.inc'
      include 'quickpot.inc'
      include 'solvents.inc'
c###############################################################
c
c     This function performs everything necessary for a solute
c     move .. returns true if the move is accepted.
c
c     (C) Christopher Woods, Julien Michel, November 2004
c
c###############################################################

      integer updatenrg,movnum,savenrg
      integer isol,ires,movtyp
      double precision ieold,ienew
      double precision soluteEnergy
      double precision extraEnergy,inrgf,inrgb
      logical metropolis
      integer lgth,lineLength
      character*4 statstring
      integer ranInt
      integer i
c     JM
      integer jsol

c     clear the move list     
      call clearMoveList()      
      
c     randomly choose a solute - all solutes are weighted evenly
      soluteMove = .false.
      if (NSolutes.le.0) return
        
100   continue
      isol = ranInt(1,NSolutes)
      if (IgnoreSolute(isol)) goto 100
c      PRINT*, "I AM",isol
c     randomly choose a residue - all residues are weighted evenly
      ires = ranInt(1,NSolResidues(isol))

c     add this solute to the move list
      call setSoluteMoved(isol)
      
c     increment the number of attempted moves for this solute
      NSolTry(isol,ires) = NSolTry(isol,ires) + 1      
            
c     see whether we are updating the solute energies
      if (updatenrg.eq.UPDATEENERGY) then
        savenrg = OLD
      else 
        savenrg = DONTCARE
      endif

c     Check if a coupled solute must be moved as well
      jsol = DualTopologySolute(isol)
      if (jsol .lt. 0) jsol = -jsol
      if (jsol .ne. 0) call setSoluteMoved(jsol)

c     get the old energy of the solute - the internal angles were calculated
c     when it was built
      ieold = soluteEnergy(savenrg,isol,inrgb,inrgf)
     .            + extraEnergy(savenrg,inrgb,inrgf)

c     JM Dual Topology get energy of old solute if it is being moved
      if (jsol .ne. 0) ieold = ieold + soluteEnergy(savenrg,jsol,inrgb,inrgf)

c     copy the coordinates to storage
      call pushBackSoluteCoords(isol)
c     JM Dual Topology consider saving coords of coupled solute
      if (jsol .ne. 0) call pushBackSoluteCoords(jsol)
      
c     move the solute
      call moveSolute(isol,ires)

c     rebuild the solute - this will recalculate the internal angles
      call buildSolute(isol)   

c     JM Dual Topology consider rebuilding coupled solute
      if (jsol .ne. 0) call buildSolute(jsol)
      
c     if we are doing FEP and have a variable solute, then synchronise
c     the forwards and backwards solutes
      if (DifferentLambdas .and. NSolVariable(isol).gt.0) then
        call synchFEPSolutes(isol)
        if (jsol .ne. 0) call synchFEPSolutes(jsol)
      endif

c     make the calculation go from old to new if necessary
      if (updatenrg.eq.UPDATEENERGY) savenrg = NEW
      
c     get the new energy
      ienew = soluteEnergy(savenrg,isol,inrgb,inrgf) + 
     .           extraEnergy(savenrg,inrgb,inrgf)

c     JM Dual Topology get energy of old solute if it is being moved
      if (jsol .ne. 0) ienew = ienew + soluteEnergy(savenrg,jsol,inrgb,inrgf)

c      write(printstring,*) 'ieold ',ieold,' ienew ',ienew
c      call printLine(DEBUG,printstring)
c     now perform the metropolis test
      if (metropolis(ieold,ienew,movtyp)) then
        soluteMove = .true.
        statstring = "PASS"            
        NSolAcp(isol,ires) = NSolAcp(isol,ires) + 1
c       update the energies
        if (updatenrg.eq.UPDATEENERGY) then
           call updateSoluteEnergies(isol)
           call updateExtraEnergies
c          JM Dual Topology, update energies of second solute if necessary
           if (jsol .ne. 0) call updateSoluteEnergies(jsol)
        endif
c       we can now clear the solute coordinate stack
        call clearSoluteCoordStack
c
c       Possible bug spotted by S. Foucher Dec. 2005
c       The prob of a solvent molecule must be updated if the solute move 
c       is accepted in case of Pref. Sampling
        if (PrefSample) then
c         we need to save the new solvent probabilities
          do i=1,NSolvents
            SvnProb(i) = SvnNewProb(i)
          enddo
          SvnMaxProb = SvnNewMaxProb
          SvnSumProb = SvnNewSumProb
        endif
c       end bugfix
      else
        soluteMove = .false.
        statstring = "FAIL"
c       copy back the old coords from the stack - this also restores
c       the internal angles
        call popSoluteCoords
c       JM this is annoying...
        if (doGB .and. .not. ConstBorn) call popPDACont
c       JM if dual Topology the stack needs to be popped twice
        if (jsol .ne. 0) call popSoluteCoords
        call clearSoluteCoordStack
      
c       if we are doing FEP and have a variable solute, then synchronise
c       the forwards and backwards solutes - this restores the old coordinates
        if (DifferentLambdas .and. NSolVariable(isol).gt.0) then
          call synchFEPSolutes(isol)
          if (jsol .ne. 0) call synchFEPSolutes(jsol)
        endif
      
      endif

c     clear the move lists if not doing a quick move
      if (.not. QuickMove) call clearMoveList()      
      
c     calculate the total energy as a check
      if (testenergy) call testEnergies
                  
c     print the status of this move
      if (printMoveEvery.eq.0) return
      if (mod(movnum,printMoveEvery).eq.0) then
        lgth = lineLength(SoluteName(isol))
        write(printstring,12) movnum,isol,SoluteName(isol)(1:lgth),SolResNam(isol,ires),
     .                        oldTotalEnergy,statstring
12      format(i10,": SOLUTEMOVE ",i5," ",a," residue ",a,1x,g20.8,1x,a4)
        call printLine(MOVE,printstring)
        if (streamhandle(ENERGY).ge.0) call printEnergies(movnum)
      endif
      
      return
      end
