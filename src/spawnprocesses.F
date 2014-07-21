      subroutine spawnProcesses(iter,doequil,doreti,dopt)
      include 'dimensions.inc'
      include 'reti.inc'
      include 'mpif.h'
c###################################################################
c
c     This routine spawns a number ProtoMS processes using
c     a slave script
c
c     Samuel Genheden, 2013
c
c###################################################################
      integer iter,doequil,doreti,dopt,ierr
      integer mpistat(MPI_STATUS_SIZE)
      integer i,j,islave,lineLength,flag
      character*300 cmd

c     Spawn new slave processes for lambdaidx > 1
      islave = 0          
      do i=2,nlambda  
        if (nptblock.eq.0) then
          write(cmd,"(2A,I6,2I3,3I2)")"python ",slave(1:lineLength(slave)),iter,i-1,0,doreti,0,doequil
          call MPI_Send(cmd,300,MPI_CHARACTER,i-1,11,MPI_COMM_WORLD,ierr)
c         write(6,*)"Sended ",cmd(1:lineLength(cmd))
        else
          do j=1,ntemp(i)
            islave = islave + 1
            flag = doreti
            if (j.gt.1) flag = 0
            write(cmd,"(2A,I6,2I3,3I2)")"python ",slave(1:lineLength(slave)),iter,i-1,j-1,flag,dopt,doequil
            call MPI_Send(cmd,300,MPI_CHARACTER,islave,11,MPI_COMM_WORLD,ierr)
          enddo
        endif
      enddo

c     Spawn new slave processes for lambdaidx = 1
      if (nptblock.gt.0) then
        do j=2,ntemp(1)
          islave = islave + 1
          write(cmd,"(2A,I6,2I3,3I2)")"python ",slave(1:lineLength(slave)),iter,0,j-1,0,dopt,doequil
          call MPI_Send(cmd,300,MPI_CHARACTER,islave,11,MPI_COMM_WORLD,ierr)
        enddo
      endif

c     Run the processes for lambdaidx = 1 and tempidx = 1 on the master node
      write(cmd,"(2A,I6,2I3,3I2)")"python ",slave(1:lineLength(slave)),iter,0,0,doreti,dopt,doequil
      call system(cmd(1:lineLength(cmd)))

c     Wait for all processes to complete
      i = 0
      do while (i.lt.nminslaves-1) 
        call MPI_Recv(cmd,300,MPI_CHARACTER,MPI_ANY_SOURCE,11,MPI_COMM_WORLD,mpistat,ierr)
c       write(6,*)"Recv ",cmd(1:lineLength(cmd))
        i = i + 1
      enddo

      end subroutine
