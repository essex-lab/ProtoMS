      subroutine splitDDDSAsimulate(nconfig,nsplit)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'pressure.inc'
      include 'flags.inc'
      include 'printer.inc'
      include 'gb.inc'
      include 'quickpot.inc'
      include 'energies.inc'
      include 'surface.inc'
c################################################################
c
c     This routine 
c     __ 1) Save coordinates and Accurate energy of a a configuration
c     __ 2) Performs split MC moves with a simple potential
c     __ 3) Accept the last generated configuration according to the 
c       energy differences between the simple and accurate potential
c     __ 4) Update observables
c     __ 5) Return to 1) or exit after nconfig repetitions
c     
c     Remarks :
c     SplitSimulate only supports residue,solute and solvent moves
c
c     Expects the correct GB potential to run with SASA on
c     (c) Julien Michel Jan 2006
c
c#################################################################

      integer nconfig,nsplit
      integer i,j,imov     
      integer NAcceptedMoves,NTotalMoves
      integer NAcceptedQuickMoves,NQuickMoves
      integer getMoveType
      logical residueMove,solventMove,soluteMove
      logical templ
      double precision multipleMovesEnergy
      double precision oldaccnrg,newaccnrg,newmovednrg,oldmovednrg
      double precision nrgb,nrgf
      double precision diffacc,diffquick
      logical accept
      logical metropolis

c     return if there are no configurations
      if (nconfig.le.0) return

      DDDon = .true.
      NQuickMoves = 0
      NAcceptedQuickMoves = 0
      NAcceptedMoves = 0
      NTotalMoves = 0
c     *** Initialisation section ***
c     get the initial state energy with correct potential
c     *** true ***
      doGB = .true.
      doSA = .true.
      call checkGotTotal
c     Save all these energies
      call StoreEnergies(ACCNRGI)
c     Now calculate approximate initial energies
      doGB = .false.
      doSA = .false.
      QuickMove = .true.
      GotInitialTotalEnergy = .false.
      call checkGotTotal
c     *** Save the Energies of Uquick(i)
      call StoreEnergies(QUICKNRGI)
c     Save the coordinates
      call StoreCoordinates(CONFI)
c     @@@@ Main Simulation @@@
      do i=1,nconfig
        NTotalMoves = NTotalMoves + 1
c       Perform nsplit MC moves
        doGB = .false.
        doSA = .false.
        QuickMove = .true.
c       Initialise split Move list
        call clearSplitMoveList
c        write(printstring,*) 'PERFORMING SPLIT MOVES DDDSA'
c        call printLine(DEBUG,printstring)
        do j=1,nsplit
c         obtain the move type based on probabilities entered by the user
          NQuickMoves = NQuickMoves + 1
          imov = getMoveType()
c          write(printstring,*) 'MOVE ',j,' imov ',imov
c          call printLine(DEBUG,printstring)
          if (imov.eq.TRYPROTEIN) then
            templ = residueMove(UPDATEENERGY,i,0)
          else if (imov.eq.TRYSOLUTE) then
            templ = soluteMove(UPDATEENERGY,i,0)
          else if (imov.eq.TRYSOLVENT) then
            templ = solventMove(UPDATEENERGY,i,0)
          else
            write(printstring,*) "Problems sorting out the move type! ",imov
            call printLine(WARNING,printstring)
            call printLine(WARNING,"This is probably a programming bug!")
            templ = .false.
          endif
c          write(printstring,*) 'Quick move ',j,' accepted ? ',templ
c          call printLine(DEBUG,printstring)
c         If the move is accepted, add what was moved to a list of 
c         moved things
          if (templ) then 
            NAcceptedQuickMoves = NAcceptedQuickMoves + 1            
            call UpdateSplitMoveList(imov)
          endif
c         NOW CLEAR THE STANDARD MOVE LIST 
          call clearMoveList()
        enddo
c       *** Save the Energies of Uquick(j)
c           in case the move is accepted
        call StoreEnergies(QUICKNRGJ)
        call StoreCoordinates(CONFJ)
c       *** Mark as Moved all the stuff that was on a SplitMoveList
        call FlagMovedBits
c       Switch to the accurate potential 
c       *** true ***
        doGB = .true.
        doSA = .true.
        QuickMove = .false.
c       *** __ revert to old coordinates i and energies ?***
        call revertCoordinates(CONFI)
c        write(printstring,*) 'BEFORE REVERT OLDSA ',OldSAEnergy
c        call printLine(DEBUG,printstring)
        call revertEnergies(ACCNRGI)
c        write(printstring,*) 'AFTER REVERT OLDSA ',OldSAEnergy
c        call printLine(DEBUG,printstring)
c       *** Evaluate the Energy of Uacc(j). 
        oldaccnrg = OldTotalEnergy
        oldmovednrg = multipleMovesEnergy(OLD,nrgb,nrgf)
c       *** __ revert to new coordinates j****
        call revertCoordinates(CONFJ)
        newmovednrg = multipleMovesEnergy(NEW,nrgb,nrgf)
        newaccnrg = oldaccnrg + newmovednrg - oldmovednrg
c        write(printstring,*) 'OLDACC ',oldaccnrg,' NEWACC ',newaccnrg,
c     .    ' NEWMOVED ',newmovednrg,' OLDMOVED ',oldmovednrg
c        call printLine(DEBUG,printstring)
c        write(printstring,*) ' NewSAEnergy ',NewSAEnergy,' F ',NewSAEnergyF,
c     .    ' B ',NewSAEnergyB
c        call printLine(DEBUG,printstring)
c        write(printstring,*) ' OldSAEnergy ',OldSAEnergy,' F ',OldSAEnergyF,
c     .    ' B ',OldSAEnergyB
c        call printLine(DEBUG,printstring)        

c        GotInitialTotalEnergy = .false.
c        call checkGotTotal

c       The correct new energies should be in place now
c        write(printstring,*) ' ENERGIES I KNOW ABOUT '
c        call prettyLine(DEBUG,printstring)
c        write(printstring,*) 'Uacc(i) ',oldaccnrg
c        call printLine(DEBUG,printstring)
c        write(printstring,*) 'Uacc(j) ',newaccnrg
c        call printLine(DEBUG,printstring)
c        write(printstring,*) 'Uquick(i) ',StoredTotalEnergy(QUICKNRGI)
c        call printLine(DEBUG,printstring)
c        write(printstring,*) 'Uquick(j) ',StoredTotalEnergy(QUICKNRGJ)
c        call printLine(DEBUG,printstring)
c       *** Perform Acceptance Test
        diffacc = newaccnrg - oldaccnrg
        diffquick = StoredTotalEnergy(QUICKNRGJ) - StoredTotalEnergy(QUICKNRGI)
        write(printstring,*) ' diffquick ',diffquick,' diffacc ',diffacc,
     .    ' delta ',diffacc-diffquick
        call printLine(DEBUG,printstring)
        accept = metropolis(diffquick,diffacc,0)
        if (accept) then 
c         New must become old...
          call updateMultipleMovesEnergies()
          write(printstring,*) '$$$$$$$SPLIT MOVE ACCEPTED$$$$$$$$'
          call printLine(DEBUG,printstring)
          NAcceptedMoves = NAcceptedMoves + 1
          call StoreEnergies(ACCNRGI)
          call StoreCoordinates(CONFI)
        else
          write(printstring,*) '$$$$$$$SPLIT MOVE REJECTED$$$$$$$'
          call printLine(DEBUG,printstring)
          call RevertEnergies(ACCNRGI)
          call RevertCoordinates(CONFI)
        endif
c       *** Clear the Moved List and residue energies stack
        call clearMoveList()
        call clearResidueEnergyStack()
c       now average the energies
c        write(printstring,*) 'Tot PROT nrg comp bnd ',OldProBnd(1),
c     .    ' ang ',OldProAng(1),' dih ',OldProDih(1),' nbc ',OldProNBC(1),
c     .    ' nbl ',OldProNBLJ(1)
c        call printLine(DEBUG,printstring)
        call calcAverageEnergy
        if (accept) then
          call RevertEnergies(QUICKNRGJ)
          call StoreEnergies(QUICKNRGI)
        else
          call popPDAcont()
          call RevertEnergies(QUICKNRGI)
        endif
      enddo
c     Disable quick move, should we be doing other things after the block 
c     of simulation (for example a RETI move)
      QuickMove = .false.
      DDDon = .false.
c     *** true ***
      doGB = .true.
      doSA = .true.
c      write(printstring,*) ' ### TESTING ENERGY CONSERVATION ### '
c      call printLine(DEBUG,printstring)
c      GotInitialTotalEnergy = .false.
c      call checkGotTotal

      if (NTotalMoves.gt.0) then
        write(printstring,*) "Data collection complete. Attempted/accepted quick moves ",NQuickMoves,NAcceptedQuickMoves
        call printLine(INFO,printstring)
        write(printstring,*) "Attempted/accepted split moves ",NTotalMoves,NAcceptedMoves
        call printLine(INFO,printstring)
        call printLine(SPLITPOT,printstring)
        
        if (streamhandle(ACCPT).ge.0) call printAttemptAccept
      endif
      
      if (SolventBox) then
        write(printstring,*) "Final volume = ",volume," A^3"
        call prettyLine(INFO,printstring)
      endif

      return
      end
