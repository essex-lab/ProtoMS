      logical function getAppendFlag(stream) 
      include 'dimensions.inc'
      include 'printer.inc'
c#####################################################
c
c     This routine sets the stream 'stream' to 
c     output to the file called filename.
c
c     (C) Christopher Woods
c
c#####################################################
      integer stream

      if (appendstreams) then
        getAppendFlag = (.not.(stream.eq.RESTART))
      else
        getAppendFlag = .false.
      endif

      return 

      end function
