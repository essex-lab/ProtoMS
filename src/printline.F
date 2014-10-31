      subroutine printLine(stream,line)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'parallel.inc'
c################################################################
c
c     This routine prints line 'line' to stream 'stream'
c
c     (C) Christopher Woods, October 2004
c
c################################################################

      integer stream
      character*(*) line
      integer lgth,lineLength,FILE,opfile

c     get the file that this stream is connected to
      FILE = streamhandle(stream)      
c     don't output the line if this stream is not connected to a file
      if (FILE.lt.0) return
      
      lgth = lineLength(line)
      if (lgth.eq.0) return
      
      if (FILE.eq.STDOUT .or. FILE.eq.STDERR) then
        write(FILE,10) streamtag(stream),line(1:lgth)
10      format(a8,a)
      else
        opfile = FILE-40
c       if the file has not been opened, then try to open it!        
        if (.not.fileopened(opfile)) then
          if (fileappend(opfile)) then
            open(unit=FILE,form="formatted",status="old",position='APPEND',err=300,file=openfiles(opfile))
          else
300         continue
            open(unit=FILE,form="formatted",status="unknown",err=100,file=openfiles(opfile))
          endif
          fileopened(opfile) = .true.
          goto 200
c-------------- This block is only executed if the file refuses to open! --------------------
100       continue
          FILE = STDOUT
          streamhandle(stream) = STDOUT
          lgth = lineLength(openfiles(opfile))
          write(FILE,10) streamtag(stream),"Error opening file ",openfiles(opfile)(1:lgth)
          write(FILE,10) streamtag(stream),"Text to this stream will be redirected to STDOUT"
          lgth = lineLength(line)
          write(FILE,10) streamtag(stream),line(1:lgth)
          return
c--------------------------------------------------------------------------------------------
200       continue
        endif

        write(FILE,11) line(1:lgth)
11      format(a)
      endif

c     flush the output (since generally the output is quite slow)
c     **THIS IS NOT A STANDARD INTRINSIC FUNCTION!**
c     ** Comment out this line if your compiler does not support it! **
      call flush(FILE)

      return
      
      end
