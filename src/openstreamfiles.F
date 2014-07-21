      subroutine openStreamFiles
      include 'dimensions.inc'
      include 'printer.inc'
      include 'splitter.inc'
c###############################################
c
c     This function opens all of the files
c     associated with the streams
c
c     (C) Christopher Woods
c
c###############################################

      integer i,j,lgth,lineLength
      logical getOption
      integer openFile
      character*14 optn
      
      do i=1,nstreams
c       write the name of the option - it will be 'stream???'
c       where ??? is the tag for the string, e.g. streaminfo or
c       streamfatal        
        call split(streamtag(i))
        lgth = lineLength(words(1))
        optn(1:6) = 'stream'
        optn(7:6+lgth) = words(1)(1:lgth)
        if (7+lgth.lt.15) then
          do j=7+lgth,14
            optn(j:j) = ' '
          enddo
        endif
        
        call lowercase(optn)
        
c       now see if this option is specified
        if (getOption(optn)) then
          streamhandle(i) = openFile(words(1))
c         ensure that the fatal and warning streams are directed to something          
          if (streamhandle(i).lt.0) then
            if (i.eq.WARNING .or. i.eq.FATAL) then
              streamhandle(i) = STDERR
            endif
          endif
        endif
      enddo
      
      return
      end
