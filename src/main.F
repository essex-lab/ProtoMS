      program main
      include 'protoms.inc'
      include 'mpif.h'
C     (C) Christopher Woods  2002-2004
c      July 2014 - Samuel Genheden added MPI stuff
c      March 2015 - Modified MPI Ana I. Cabedo Martinez

      integer ierr

      mylamid = 0
      mytemid = 0
      mygcid = 0
      myjawsid = 0

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,ntasks,ierr)

c     initialise the print system
      call initialisePrintSystem

c     parse any command line arguments
      call parseArgs

c     setup replica exchange
      call setupReplicas

      write(6,*)myid," MAIN SUR mt_beta_t",mt_beta_t

c     set the default locations of any streams
      call openStreamFiles

c     initialise the program
      call initialise
      
c     set the system up
      call setup

      write(6,*)myid," MAIN GO mt_beta_t",mt_beta_t
      write(6,*)myid," MAIN mylamid ",mylamid," mytemid ",mytemid
      
c     run the simulation
      call simulation
      
c     shut everything down / print out restart/info 
      call shutdown

      call MPI_FINALIZE(ierr)

      end
