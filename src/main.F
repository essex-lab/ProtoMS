      program main
      include 'protoms.inc'
      include 'mpif.h'
C     (C) Christopher Woods  2002-2004
c      July 2014 - Samuel Genheden added MPI stuff

      integer ierr

      mylamid = mytemid = mygcid = myjawsid = 0
      nlamtasks = ntemtasks = ngctasks njawstasks = 1

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,ntasks,ierr)

c     initialise the print system
      call initialisePrintSystem

c     parse any command line arguments
      call parseArgs

c     setup replica exchange
      call setupReplicas

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

      call MPI_FINALIZE(ierr)

      end
