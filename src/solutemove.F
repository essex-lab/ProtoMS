      
      logical function soluteMove(updatenrg,movnum,movtyp)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'printer.inc'
      include 'solutes.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'moves.inc'
      include 'energies.inc'
      include 'gb.inc'
      include 'quickpot.inc'
      include 'solvents.inc'
      include 'proteins.inc'
      include 'simparams.inc'
c###############################################################
c
c     This function performs everything necessary for a solute
c     move .. returns true if the move is accepted.
c
c     (C) Christopher Woods, Julien Michel, November 2004
c     Modified by Samuel Genheden, 2014
c
c###############################################################

      integer updatenrg,movnum,savenrg
      integer isol,ires,movtyp
      double precision ieold,ienew,penalty
      double precision soluteEnergy
      double precision extraEnergy,inrgf,inrgb
      logical metropolis
      integer lgth,lineLength
      character*4 statstring
      integer ranInt
      integer i, j
c     JM
      integer jsol, idih

      double precision tunesolute, dihmovesize
      integer istrt,iend
      
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

c     In group dual topology each solute is moved independently
      if (NSolGroups.gt.0) jsol = 0

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
      
      if (tunedihedral) then
c     pick a flexible dihedral  at random
         iDih = ranInt ( 1, NRotDihs(isol,ires) )
c     move the solute
         dihmovesize = tuneSolute(isol,ires,iDih)
      else
         call moveSolute(isol,ires)
         iDih = 0
         dihmovesize = 0
      endif

c     if this is theta solute, it has to reside within a box
c     check that solute is inside the given box
      penalty = ZERO
      if ((NThetaSolutes.gt.0).and.(isol.le.NThetaSolutes)) then
        if (SolGeomCentre(isol,1).gt.(coordx+gridx+TINY)) then
          penalty = (HUGE*HUGE)
        endif
        if (SolGeomCentre(isol,2).gt.(coordy+gridy+TINY)) then
          penalty = (HUGE*HUGE)
        endif
        if (SolGeomCentre(isol,3).gt.(coordz+gridz+TINY)) then
          penalty = (HUGE*HUGE)
        endif
        if (SolGeomCentre(isol,1).lt.(coordx-TINY)) then
           penalty = (HUGE*HUGE)
        endif
        if (SolGeomCentre(isol,2).lt.(coordy-TINY)) then
           penalty = (HUGE*HUGE)
        endif
        if (SolGeomCentre(isol,3).lt.(coordz-TINY)) then
          penalty = (HUGE*HUGE)
        endif
      endif

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
     .           extraEnergy(savenrg,inrgb,inrgf) + penalty

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
c       The prob of all solvent and protein molecules must be updated if the solute move 
c       is accepted in case of Pref. Sampling
        if (PrefSample) then
           do i=1,NSolvents
              SvnProb(i) = SvnNewProb(i)
           enddo
           SvnMaxProb = SvnNewMaxProb
           SvnSumProb = SvnNewSumProb
        endif
        if (NProteins.gt.0 .and. ProPrefSample) then
           do i=1,NProteins
              do j=1,NResidues(i)
                 ProProb(i,j) = ProNewProb(i,j)
              enddo
           enddo
           ProMaxProb = ProNewMaxProb
           ProSumProb = ProNewSumProb
        endif
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

c     if we are doing a dihedral tuning simulation then record statistics of
c     move successes for individual dihedrals
      if (tunedihedral) then
         do i=1,20
            if (abs(dihmovesize)/PI*180.0.lt.float(i)) then
               DihMoveAttempts(isol,ires,iDih,i) = DihMoveAttempts(isol,ires,iDih,i) + 1
               if (soluteMove) DihMoveSuccess(isol,ires,iDih,i) = DihMoveSuccess(isol,ires,iDih,i) + 1
               goto 200
            endif
         enddo
      endif
      
c     clear the move lists if not doing a quick move
 200  if (.not. QuickMove) call clearMoveList()      
      
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
