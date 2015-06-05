      program main
      include 'protoms.inc'
#ifdef WITH_MPI_
      include 'mpif.h'
#endif
C     (C) Christopher Woods  2002-2004
c      July 2014 - Samuel Genheden added MPI stuff
c      March 2015 - Modified MPI Ana I. Cabedo Martinez

      integer ierr

      mylamid = 0
      mytemid = 0
      mygcid = 0
      myjawsid = 0

#ifdef WITH_MPI_
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,ntasks,ierr)
#endif

c     initialise the print system
      call initialisePrintSystem

c     parse any command line arguments
      call parseArgs

c     setup replica exchange
#ifdef WITH_MPI_
      call setupReplicas
#endif

c     set the default locations of any streams
      call openStreamFiles

c     initialise the program
      call initialise
      
c     set the system up
      call setup
      
c     run the simulation
      call simulation
      
c     shut everything down / print out restart/info 
      call shutdown

#ifdef WITH_MPI_
      call MPI_FINALIZE(ierr)
#endif

      end
