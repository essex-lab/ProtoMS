      subroutine closeProgram(stat,str)
      implicit none
      include 'printer.inc'
#ifdef WITH_MPI_
      include 'mpif.h'
#endif
c#######################################################
c
c     This closes the program because of an error
c     May eventually drop different error codes to
c     signify different behaviour...
c
c     (C) Christopher Woods    30 January 2002
c
c#######################################################

      character*(*) str
      integer stat,ierr
      
      call printLine(FATAL,str)
      write(printstring,*) "Program exiting with status ",stat      
      call printLine(FATAL,printstring)
            
c     try to print out results and restart information
      call shutdown()
      flush(FATAL)

#ifdef WITH_MPI_
      call MPI_FINALIZE(ierr)
#endif
      call exit(stat)
      
      end
