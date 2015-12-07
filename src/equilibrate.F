      subroutine equilibrate(nconfig)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
c###########################################################
c
c     This routine performs the equilibration of the 
c     entire system
c
c     (C) Christopher Woods     9 August 2002
c
c###########################################################
                                     
      double precision tstart,tend,telapsed
      integer nconfig
      integer i,imov     
      integer NAcceptedMoves,NTotalMoves
      integer getMoveType
      logical residueMove,titrationMove,solventMove,soluteMove
      logical thetaMove, sampleMove,gcsoluteMove,deletionMove,insertionMove
      logical volumeMove,rotamerMove,solThetaMove
      logical templ
                     
c     return if there are no configurations
      if (nconfig.le.0) return

c     check that the initial total energy has been calculated
      call checkGotTotal
      
      NAcceptedMoves = 0
      NTotalMoves = 0

      call cpu_time(tstart)
      
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
          write(printstring,*) "Problems sorting out the move type in equilibrate()!"
          call closeProgram(-1,printstring)
          templ = .false.
        endif
        if (templ) NAcceptedMoves = NAcceptedMoves + 1

      enddo
            
      call cpu_time(tend)
      telapsed = tend - tstart
    
      write(printstring,*) "Equilibration complete. Attempted/accepted moves ",
     .                           NTotalMoves,NAcceptedMoves
      call printLine(INFO,printstring)

      write(printstring,"('These moves took ',EN8.1E1,' s to run.')") telapsed
      call printLine(INFO,printstring)

      return
      end
