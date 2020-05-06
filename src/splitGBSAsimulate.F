      subroutine splitGBSAsimulate(nconfig,nsplit)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'pressure.inc'
      include 'flags.inc'
      include 'printer.inc'
c      include 'gb.inc'
c      include 'surface.inc'
      include 'quickpot.inc'
      include 'energies.inc'
      include 'simparams.inc'
      include 'extraenergies.inc'
c################################################################
c
c     This routine 
c     is a modified version of splitsimulate and it uses 
c     a GB potential as a quick potential instead of a
c     Distance Dependent Dielectric potential
c     The differences with the exact GB potential are 
c     --> The Born radii are updated at each step and 
c     consistently so that detailed balance is achieved
c     --> A different cutoff
c     --> No calculation of forwards/backwards energies
c
c     In addition, it expects the full potential to do SA calculations 
c     (although they are not performed in the quick potential )
c
c     (c) Julien Michel January 2006
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
      double precision accut,accut2,accborncut,quickborncut
      double precision accborncut2,quickborncut2
      double precision accthreshold,quickthreshold
      double precision quickcut,quickcut2
      double precision acclam,acclamF,acclamB,quicklam,quicklamF,quicklamB

c      double precision temp,tempf,tempb
c      double precision GBEnergy
c     return if there are no configurations
      if (nconfig.le.0) return
      NQuickMoves = 0
      NAcceptedQuickMoves = 0
      NAcceptedMoves = 0
      NTotalMoves = 0
c     *** Initialisation section ***
c     get the initial state energy with correct potential
      accut = cut
      accut2 = cut2
      accborncut = BornCut
      accborncut2 = BornCut2
      accthreshold = gbthreshold
      acclam = Lambda
      acclamF = LambdaF
      acclamB = LambdaB
      doSA = .true. 
c     Quick sim parameters
      quickborncut = 12.0d0
      quickborncut2 = quickborncut**2
      quickcut = 6.0d0
      quickcut2 = quickcut**2
      quickthreshold = 0.05
      quicklam = Lambda
      quicklamF = Lambda
      quicklamB = Lambda
C     *** Active for DEBUGGING ***
c$$$      quickborncut = accborncut
c$$$      quickborncut2 = accborncut2
c$$$      quickcut = accut
c$$$      quickcut2 = accut2
c$$$      quickthreshold = accthreshold
      call checkGotTotal
c      write(printstring,*) '1 CALLED TOTAL'
c      call printLine(WARNING,printstring)
c     Save all these energies
      call StoreEnergies(ACCNRGI)
      call StoreBornData(ACCNRGI)
c     Now calculate approximate initial energies
      QuickMove = .true.
      doSA = .false.
      cut = quickcut
      cut2 = quickcut2
      BornCut = quickborncut
      BornCut2 = quickborncut2
      gbthreshold = quickthreshold
      Lambda = quicklam
      LambdaF = quicklamF
      LambdaB = quicklamB
      DifferentLambdas = .false.
      GotInitialTotalEnergy = .false.
      call checkGotTotal
c      write(printstring,*) '1 CALLED TOTAL'
c      call printLine(WARNING,printstring)

c     *** Save the Energies of Uquick(i)
      call StoreEnergies(QUICKNRGI)
      call StoreBornData(QUICKNRGI)
c     Save the coordinates
      call StoreCoordinates(CONFI)
c     @@@@ Main Simulation @@@
      do i=1,nconfig
        NTotalMoves = NTotalMoves + 1
c       Perform nsplit MC moves
        QuickMove = .true.
        doSA = .false.
        cut = quickcut
        cut2 = quickcut2
        BornCut = quickborncut
        BornCut2 = quickborncut2
        gbthreshold = quickthreshold
        Lambda = quicklam
        LambdaF = quicklamF
        LambdaB = quicklamB
        DifferentLambdas = .false.
c       Initialise split Move list
        call clearSplitMoveList
c        write(printstring,*) 'PERFORMING SPLIT MOVES GBSA'
c        call printLine(DEBUG,printstring)
c        write(printstring,*) ' OldSAEnergy ',OldSAEnergy
c        call printLine(DEBUG,printstring) 
c        if (i .eq. 2) then 
c          write(printstring,*) ' ### TESTING ENERGY CONSERVATION ### '
c          call printLine(DEBUG,printstring)
c          GotInitialTotalEnergy = .false.
c          call checkGotTotal            
c          write(printstring,*) ' GBEnergy should be ',OldGBEnergy
c          call printLine(DEBUG,printstring) 
c        endif
        do j=1,nsplit
c          write(printstring,*) 'Born 20',ZSolEffBorn(1,20),' Old ',OldZSolEffBorn(1,20)
c          call printLine(DEBUG,printstring)
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
            write(printstring,"(a,i4)") "Problems sorting out the move type! ",imov
            call printLine(WARNING,printstring)
            call printLine(WARNING,"This is probably a programming bug!")
            templ = .false.
          endif
c          write(printstring,*) 'Born 20',ZSolEffBorn(1,20),' Old ',OldZSolEffBorn(1,20)
c          call printLine(DEBUG,printstring)
c          write(printstring,*) 'Quick move ',j,' Nrg ',OldTotalEnergy,
c     .      ' accepted ? ',templ
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
        call StoreBornData(QUICKNRGJ)
c       *** Mark as Moved all the stuff that was on a SplitMoveList
        call FlagMovedBits
c       Switch to the accurate potential 
c       *** true ***
        QuickMove = .false.
        doSA = .true.
        cut = accut
        cut2 = accut2
        BornCut = accborncut
        BornCut2 = accborncut2
        gbthreshold = accthreshold
        Lambda  = acclam 
        LambdaF = acclamF
        LambdaB = acclamB 
        DifferentLambdas = .true.
c       *** __ revert to old coordinates i and energies ?***
        call revertCoordinates(CONFI)
        call revertBornData(ACCNRGI)
        call revertEnergies(ACCNRGI)
c        write(printstring,*) 'REVERTED OLDSA ',OldSAEnergy, ' XTRA ',OldExtraEnergy
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

c       The correct new energies should be in place now
c        write(printstring,*) ' ENERGIES I KNOW ABOUT '
c        call prettyLine(DEBUG,printstring)
c        write(printstring,*) 'Uacc(i) ',oldaccnrg,' GB ',OldGBEnergy,
c     .    ' Xtra ',oldextraEnergyPart
c        call printLine(DEBUG,printstring)
c        write(printstring,*) 'Uacc(j) ',newaccnrg, ' GB ',NewGBEnergy,
c     .    ' Xtra ',newextraEnergyPart
c        call printLine(DEBUG,printstring)
c        write(printstring,*) 'Uquick(i) ',StoredTotalEnergy(QUICKNRGI),
c     .    ' GB ',StoredGBEnergy(QUICKNRGI),' Xtra ',StoredExtraEnergy(QUICKNRGI)
c        call printLine(DEBUG,printstring)
c        write(printstring,*) 'Uquick(j) ',StoredTotalEnergy(QUICKNRGJ),
c     .    ' GB ',StoredGBEnergy(QUICKNRGJ),' Xtra ',StoredExtraEnergy(QUICKNRGJ)
c        call printLine(DEBUG,printstring)

c       *** Perform Acceptance Test
        diffacc = newaccnrg - oldaccnrg
        diffquick = StoredTotalEnergy(QUICKNRGJ) - StoredTotalEnergy(QUICKNRGI)
c        write(printstring,*) ' diffquick ',diffquick,' diffacc ',diffacc,
c     .    ' delta ',diffacc-diffquick
c        call printLine(DEBUG,printstring)
        accept = metropolis(diffquick,diffacc,0)
        if (accept) then 
c         New must become old...
          call updateMultipleMovesEnergies()
c          write(printstring,*) '$$$$$$$SPLIT MOVE ACCEPTED$$$$$$$$',newaccnrg,NewSAEnergy
c          call printLine(DEBUG,printstring)
          NAcceptedMoves = NAcceptedMoves + 1
          call StoreEnergies(ACCNRGI)
          call StoreCoordinates(CONFI)
          call StoreBornData(ACCNRGI)
        else
c          write(printstring,*) '$$$$$$$SPLIT MOVE REJECTED$$$$$$$',oldaccnrg,OldSAEnergy
c          call printLine(DEBUG,printstring)
          call RevertEnergies(ACCNRGI)
          call RevertCoordinates(CONFI)
          call RevertBornData(ACCNRGI)
        endif
c       *** Clear the Moved List and residue energies stack
        call clearMoveList()
        call clearResidueEnergyStack()
c       now average the energies
        call calcAverageEnergy
        if (accept) then
          call RevertEnergies(QUICKNRGJ)
          call RevertBornData(QUICKNRGJ)
cc         The GB energy is wrong because it matches the constant born radii used before 
cc         they were updated. So I need to calculate a new gb energy using the quickGB
cc         and save Old=New
c          QuickMove = .true.
c          cut = quickcut
c          cut2 = quickcut2
c          DifferentLambdas = .false.
c          call SwitchBornRadii(CONFI,CONFJ)
cc          write(printstring,*) 'Born 20',ZSolEffBorn(1,20),' Old ',OldZSolEffBorn(1,20)
cc          call printLine(DEBUG,printstring)
c          temp = GBEnergy(OLD, tempf, tempb)
c          temp = GBEnergy(NEW, tempf, tempb)
c          write(printstring,*) 'GB corrected from ',OldGBEnergy,' to ',NewGBEnergy
c          call printLine(DEBUG,printstring)
c          OldGBEnergy = NewGBEnergy
c          OldGBEnergyF = NewGBEnergyF
c          OldGBEnergyB = NewGBEnergyB
c          call updateBornRadii()
          call StoreEnergies(QUICKNRGI)
        else
          call popPDAcont()
          call RevertEnergies(QUICKNRGI)
          call RevertBornData(QUICKNRGI)
        endif
      enddo
c     Disable quick move, should we be doing other things after the block 
c     of simulation (for example a RETI move)
      QuickMove = .false.
      doSA = .true. 
      cut = accut
      cut2 = accut2
      BornCut = accborncut
      BornCut2 = accborncut2
      gbthreshold = accthreshold
      Lambda = accLam
      LambdaF = accLamF
      LambdaB = accLamB
      DifferentLambdas = .true.
      
c      write(printstring,*) ' ### TESTING ENERGY CONSERVATION ### '
c      call printLine(DEBUG,printstring)
c      GotInitialTotalEnergy = .false.
c      call checkGotTotal

      if (NTotalMoves.gt.0) then
        write(printstring,"(a,2i12)")
     .        "Data collection complete. Attempted/accepted quick moves ",NQuickMoves,NAcceptedQuickMoves
        call printLine(INFO,printstring)
        write(printstring,"(a,2i12)") "Attempted/accepted split moves ",NTotalMoves,NAcceptedMoves
        call printLine(INFO,printstring)
        call printLine(SPLITPOT,printstring)
        
        if (streamhandle(ACCPT).ge.0) call printAttemptAccept
      endif
      
      if (SolventBox) then
        write(printstring,"(a,f12.4,a)") "Final volume = ",volume," A^3"
        call prettyLine(INFO,printstring)
      endif

      return
      end
