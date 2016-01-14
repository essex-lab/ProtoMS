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
#ifdef WITH_MPI_
      include 'mpif.h'
#endif
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
      integer i,imov,j,k,ierr
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
          fdumps(ndumps) = nint(getDouble(words(1),-1))
          if (fdumps(ndumps).gt.0) then
            kdumps(ndumps) = k
          else
            ndumps = ndumps -1
            write(printstring,"(2a)") "Ignoring dump: unable to parse frequency or frequency less than 1: ", words(1)
            call printLine(WARNING,printstring)            
          endif  
        endif
      enddo

      append = .TRUE.

      NAcceptedMoves = 0
      NTotalMoves = 0
      do i=1,nconfig
        if (myid.eq.1) write(6,*)myid,"###########################################",i
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
          if (myid.eq.1) write(6,*) "Doing insertion move."
          templ = insertionMove(UPDATEENERGY,i,0)
        else if (imov.eq.TRYDELETION) then
          if (myid.eq.1) write(6,*) "Doing deletion move."
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
          write(printstring,"(a,i4)") "Problems sorting out the move type! ",imov
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

#ifdef WITH_MPI_
c       Do lambda replica exchange
        if (doLambdaRE.and.(lrefreq.gt.0).and.(i.gt.1).and.(mod(i,lrefreq).eq.0)) then
          call lambdaRE
          if (doTemperatureRE) then
            call assigntemperaturere(mytemid,.true.,.false.)
          else if (doMultiGCMC) then
            call assignadamre(mygcid,.true.,.false.)
          endif
        endif

c       Do temperature replica exchange
        if (doTemperatureRE.and.(trefreq.gt.0).and.(i.gt.1).and.(mod(i,trefreq).eq.0)) then
c          write(6,*)myid, " calling TRE. Move ",i
          call temperatureRE
          if (doLambdaRE) then
            call assignlambdare(mylamid,.true.,.false.)
          else if (doMultiGCMC) then
            call assignadamre(mygcid,.true.,.false.)
          endif
        endif

c		Do GCMC replica exchange
        if (doMultiGCMC.and.(gcrefreq.gt.0).and.(i.gt.1).and.(mod(i,gcrefreq).eq.0)) then
c       GR This if statement is will help track which nodes are going where
          call gcmcRE
          if (doLambdaRE) then
            call assignlambdare(mylamid,.true.,.false.)
          else if (doTemperatureRE) then
            call assigntemperaturere(mytemid,.true.,.false.)
          endif 
        endif
#endif
        
      enddo

      if (NTotalMoves.gt.0) then
        write(printstring,"(a,2i12)") "Data collection complete. Attempted/accepted moves ",
     .                             NTotalMoves,NAcceptedMoves
        call printLine(INFO,printstring)
        
        if (streamhandle(ACCPT).ge.0) call printAttemptAccept
      endif

      if (doLambdaRE.and.(lrefreq.gt.0).and.mylamladder.ne.-1) then       
        write(printstring,"(A,I3,A,2I6,A,F8.3,A)")
     .     "Attempted/successful lambda swaps for lambdaladder ",mylamladder," : ",
     .     nlretry,nlreswap," (",dble(nlreswap)/dble(nlretry)*100,"%)" 
        call printLine(INFO,printstring)
      endif

      if (doTemperatureRE.and.(trefreq.gt.0).and.mytemladder.ne.-1) then       
        write(printstring,"(A,I3,A,2I6,A,F8.3,A)")
     .     "Attempted/successful temperature swaps for temperatureladder ",mytemladder," : ",
     .     ntretry,ntreswap," (",dble(ntreswap)/dble(ntretry)*100,"%)" 
        call printLine(INFO,printstring)
      endif

      if (doMultiGCMC.and.(gcrefreq.gt.0)) then   
c        write(6,*)"### GR simulate.F: Printing acceptance"    
        write(printstring,"(A,I3,A,2I6,A,F8.3,A)")
     .     "Attempted/successful GCMC swaps for ID ",myid," : ",
     .     ngcretry,ngcreswap," (",dble(ngcreswap)/dble(ngcretry)*100,"%)" 
        call printLine(INFO,printstring)
      endif
     
      if (SolventBox) then
        write(printstring,"(a,f12.4,a)") "Final volume = ",volume," A^3"
        call prettyLine(INFO,printstring)
      endif

      return
      end
