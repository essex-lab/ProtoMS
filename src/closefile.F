      subroutine closeFile(filename)
      include 'dimensions.inc'
      include 'printer.inc'
c################################################
c
c     This routine closes the file called
c     filename if it is currently open
c
c     (C) Christopher Woods, October 2004
c
c################################################

      character*(*) filename
      integer lgth,lineLength
      integer i,file
      
      lgth = lineLength(filename)
      if (lgth.le.0) return
      
c     see if this file has already been opened
      do i=1,nOpenFiles
        if (filename(1:lgth).eq.openfiles(i)(1:lgth)) then
c         the file is assigned to a unit
          file = i+40
          if (fileopened(i)) then
c           the file is opened - close it!
            close(file)
            fileopened(i) = .false.
          endif          
          return
        endif
      enddo
            
      return
      end
