      subroutine killSlaves(ntasks)
      include 'dimensions.inc'
      include 'reti.inc'
      include 'mpif.h'
      integer ntasks,i,ierr
      character*300 cmd

      cmd = "Exit"
      do i=1,ntasks-1  
        call MPI_Send(cmd,300,MPI_CHARACTER,i,11,MPI_COMM_WORLD,ierr)
      enddo

      end subroutine
