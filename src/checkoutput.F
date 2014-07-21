      subroutine checkOutput(stat)

      include 'dimensions.inc'
      include 'reti.inc'
      include 'mpif.h'
      include 'printer.inc'
c###################################################################
c
c     This routine check that protoMS has created all restart files
c
c     Samuel Genheden, 2013
c
c###################################################################
      integer i,j,w,stat,test,lineLength
      logical fileExists

      w = 0 
      stat = 0
      do while (w.lt.2000) 
        test = 0
        do i=1,nlambda
          if (nptblock.eq.0) then
            if (fileExists(rstfiles(i,1))) test = test + 1               
          else
            do j=1,ntemp(i)
              if (fileExists(rstfiles(i,j))) test = test + 1
            enddo
          endif        
        enddo     
        if (test.lt.nminslaves) then
          call system("sleep 4s")   
          w = w + 1
        else
          w = 4000
        endif
      enddo

      if (test.lt.nminslaves) then
        write(printstring,"(A,F5.1,A)")"It seems that not all necessary files were created by ProtoMS! Exit."
        call prettyLine(FATAL,printstring) 
        do i=1,nlambda
          if (nptblock.eq.0) then
            if (.not.fileExists(rstfiles(i,1))) write(6,*)rstfiles(i,1)(1:lineLength(rstfiles(i,1)))
          else
            do j=1,ntemp(i)
              if (.not.fileExists(rstfiles(i,j))) write(6,*)rstfiles(i,j)(1:lineLength(rstfiles(i,j)))
            enddo
          endif
        enddo 
        stat = 1000
      endif

      end subroutine
