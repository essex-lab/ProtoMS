      subroutine checkSameStart
      include 'dimensions.inc'
      include 'printer.inc'
      include 'parallel.inc'
      include 'moves.inc'
      include 'mpif.h'
      
c#######################################################
c
c     This routine ensures that different processes are
c     appropriately in sync when restarting simulations.
c
c     (C) Chris Cave-Ayland    21 November 2017
c
c#######################################################

      integer starts(ntasks)
      logical completes(ntasks)
      integer ierr, i

c     gather startstep and completed values from all processes to make sure
c     everyone is starting in the same place.
      call MPI_AllGather(completed, 1, MPI_LOGICAL, completes, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)      
      call MPI_AllGather(startstep, 1, MPI_INTEGER, starts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

c     check that processes are consistently restarting a simulation or not
      do i=2,ntasks
         if (completes(i).neqv.completes(1)) then
c     raise all hell if they don't agree
            write(printstring,*) "Not all processes attempting to restart simulation. This process:",
     .           .not.completed
            call closeProgram(-1, printstring)
         endif
      enddo

c     all start/restarting together, but are all processes in the right place?
c     check all start values against one another
      do i=2,ntasks
         if (starts(i).ne.starts(1)) then
c     raise all hell if they don't agree
            write(printstring,*) "Different processes attempting to start in different places.",
     .           " This process starts at", startstep
            call closeProgram(-1, printstring)
         endif
      enddo
      end
