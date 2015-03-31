      subroutine simulate(nconfig)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'pressure.inc'
      include 'flags.inc'
      include 'printer.inc'
      include 'moves.inc'
      include 'stacks.inc'
      include 'keyset.inc'
      include 'splitter.inc'
      include 'parallel.inc'
      include 'averages.inc'
      include 'mpif.h'
c################################################################
c
c       This routine actually runs a simulation with some
c       data collection! It produces the averages that are
c       necessary to generate thermodynamic quantities,
c       or records energies/lambdas that are necessary for
c       umbrella type or generalised ensemble type simulations!
c
c#################################################################

      integer nconfig
      integer i,imov,j,k    
      integer NAcceptedMoves,NTotalMoves,linelength
      integer getMoveType,numb
      logical residueMove,titrationMove,solventMove,soluteMove
      logical rotamerMove
      logical volumeMove,sampleMove
      logical insertionMove,deletionMove
      logical gcsoluteMove,thetaMove,solThetaMove
      logical templ,append
      character*64 tempstring
      integer ndumps, kdumps(MAXDUMP),fdumps(MAXDUMP)
      double precision getDouble
                     
c     return if there are no configurations
      if (nconfig.le.0) return
                               
c     check that the initial total energy has been calculated
      call checkGotTotal

      ndumps = 0
      do k=1,nkeys
c       4 in the line below stands for the lenght of 'dump'
        if (keys(k)(1:4).eq.'dump') then
          ndumps = ndumps + 1
          call split(values(k))
          call lowercase(words(1))
          fdumps(ndumps) = nint(getDouble(words(1),0))
          kdumps(ndumps) = k
        endif
      enddo

      append = .TRUE.

      NAcceptedMoves = 0
      NTotalMoves = 0
      do i=1,nconfig
        NTotalMoves = NTotalMoves + 1

c       obtain the move type based on probabilities entered by the user
        imov = getMoveType()
        if (imov.eq.TRYPROTEIN) then
          templ = residueMove(UPDATEENERGY,i,0)
        else if (imov.eq.TRYSOLUTE) then
          templ = soluteMove(UPDATEENERGY,i,0)
        else if (imov.eq.TRYTITRATE) then
          templ = titrationMove(UPDATEENERGY,i,0)
        else if (imov.eq.TRYSOLVENT) then
          templ = solventMove(UPDATEENERGY,i,0)
        else if (imov.eq.TRYVOLUME) then
          templ = volumeMove(UPDATEENERGY,i,0)
        else if (imov.eq.TRYROTAMER) then
          templ = rotamerMove(UPDATEENERGY,i,0)
        else if (imov.eq.TRYINSERTION) then
          templ = insertionMove(UPDATEENERGY,i,0)
        else if (imov.eq.TRYDELETION) then
          templ = deletionMove(UPDATEENERGY,i,0)
        else if (imov.eq.TRYGCSOLUTE) then
          templ = gcsoluteMove(UPDATEENERGY,i,0)
        else if (imov.eq.TRYSAMPLE) then
          templ = sampleMove(UPDATEENERGY,i,0)
        else if (imov.eq.TRYTHETA) then
          templ = thetaMove(UPDATEENERGY,i,0)
        else if (imov.eq.TRYSOLTHETA) then
          templ = solThetaMove(UPDATEENERGY,i,0)
        else
          write(printstring,*) "Problems sorting out the move type! ",imov
          call printLine(WARNING,printstring)
          call printLine(WARNING,"This is probably a programming bug!")
          templ = .false.
        endif
        if (templ) NAcceptedMoves = NAcceptedMoves + 1

c       now average the energies
        call calcAverageEnergy

c       This we force re-calculation of BAR, MBAR and TI energies
        if (doLambdaRE) reCalculateFreeEnergies = .true.

c       taking care of dumping with the right frequency
        do j=1,ndumps
          if (mod(i,fdumps(j)).eq.0) then
            append = .not.(i.eq.fdumps(j))
            call split(values(kdumps(j)))
            call join(2,nwords,tempstring)
            numb = int(i/fdumps(j))
            call dump(tempstring,append,numb)
          endif
        enddo

c       Do lambda replica exchange
        if (doLambdaRE.and.(lrefreq.gt.0).and.(i.gt.1).and.(mod(i,lrefreq).eq.0)) then
          call lambdaRE
        endif

c       Do lambda replica exchange
        if (doTemperatureRE.and.(trefreq.gt.0).and.(i.gt.1).and.(mod(i,trefreq).eq.0)) then
          call temperatureRE
        endif
        
      enddo

      if (NTotalMoves.gt.0) then
        write(printstring,*) "Data collection complete. Attempted/accepted moves ",
     .                             NTotalMoves,NAcceptedMoves
        call printLine(INFO,printstring)
        
        if (streamhandle(ACCPT).ge.0) call printAttemptAccept
      endif

      if (doLambdaRE.and.(lrefreq.gt.0)) then       
        write(printstring,"(A,2I6,A,F8.3,A)")
     .     "Attempted/successful lambda swaps: ",nlretry,nlreswap," (",dble(nlreswap)/dble(nlretry)*100,"%)" 
        call printLine(INFO,printstring)
      endif

      if (doTemperatureRE.and.(trefreq.gt.0)) then       
        write(printstring,"(A,2I6,A,F8.3,A)")
     .     "Attempted/successful temperature swaps: ",ntretry,ntreswap," (",dble(ntreswap)/dble(ntretry)*100,"%)" 
        call printLine(INFO,printstring)
      endif
     
      if (SolventBox) then
        write(printstring,*) "Final volume = ",volume," A^3"
        call prettyLine(INFO,printstring)
      endif

      return
      end
