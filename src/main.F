      program main
      include 'protoms.inc'
      include 'mpif.h'
C     (C) Christopher Woods  2002-2004
c      July 2014 - Samuel Genheden added MPI stuff
c      March 2015 - Modified Ana I. Cabedo Martinez

      integer ierr

      mylamid = 0
      mytemid = 0
      mygcid = 0
      myjawsid = 0
      nlamtasks = 1
      ntemtasks = 1
      ngctasks = 1
      njawstasks = 1

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,ntasks,ierr)

c     initialise the print system
      call initialisePrintSystem

c     parse any command line arguments
      call parseArgs

      write(6,*)"I'm before setupreplicas"

c     setup replica exchange
      call setupReplicas

      write(6,*)"Main_Before_OpenStreamFiles ",temperatureC

c     set the default locations of any streams
      call openStreamFiles

      write(6,*)"Main_After_OpenStreamFiles ",temperatureC
      
c     initialise the program
      call initialise

      write(6,*)"Main_After_Initialize ",temperatureC
      
c     set the system up
      call setup

      write(6,*)"Main_After_SetUp ",temperatureC
      
c     run the simulation
      call simulation

      write(6,*)"Main_After_Simulation ",temperatureC
      
c     shut everything down / print out restart/info 
      call shutdown

      call MPI_FINALIZE(ierr)

      end
