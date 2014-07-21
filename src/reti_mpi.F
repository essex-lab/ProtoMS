      program RETIMPI
      include 'dimensions.inc'
      include 'reti.inc'
      include 'printer.inc'
      include 'mpif.h'
c#######################################################
c
c     Main program for RETI simulations using 
c     Open MPI
c
c     Samuel Genheden, February 2013
c
c#######################################################
       
      integer ierr,myid,ntasks

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,ntasks,ierr)

      if (myid.eq.0) then
c       initialise the print system
        call initialisePrintSystem

c       parse any command line arguments
        call parseRetiArgs

c       Write some info
        write(printstring,"(A,I5,A)")" Running on ",ntasks," processors"
        call prettyline(INFO,printstring)

c       load input options
        call getretioptions

c       do the simulations
        call retiSimulation_mpi(myid,ntasks)

      else
        call retiSimulation_mpi(myid,ntasks)
      endif

      call MPI_FINALIZE(ierr)

      end program
