      subroutine closeFiles
      include 'dimensions.inc'
      include 'printer.inc'
c################################################
c
c     This routine closes all open files
c
c     Samuel Genheden, July 2013
c
c################################################

      integer i,file
      
      
      do i=1,nOpenFiles
c       the file is assigned to a unit
        file = i+40
        if (fileopened(i)) then
c         the file is opened - close it!
          close(file)
          fileopened(i) = .false.
        endif                  
      enddo
            
      return
      end
