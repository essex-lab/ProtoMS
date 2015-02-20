      subroutine setStream(stream,filename)
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

      character*(*) stream
      character*(*) filename
      character*20 tmpstr1,tmpstr2
      integer i,openFILE,lgth,lgth2,lineLength
      logical getAppendFlag
      
      lgth = lineLength(stream)
      if (lgth.gt.20) lgth = 20
      
      tmpstr1 = stream(1:lgth)
      call lowercase(tmpstr1)
      
      do i=1,nstreams
        tmpstr2 = streamtag(i)
        call leftJustify(tmpstr2)
        call lowercase(tmpstr2)
        
        lgth2 = lineLength(tmpstr2)
        
        if (tmpstr1(1:lgth).eq.tmpstr2(1:lgth2)) then
          streamhandle(i) = openFILE(filename)
c         ensure that the FATAL and WARNING streams are always on
          if (streamhandle(i).lt.0) then
            if (i.eq.FATAL .or. i.eq.WARNING) then
              streamhandle(i) = STDERR
            endif
          else
            fileappend(streamhandle(i)-40) = getAppendFlag(i)
          endif
          return
        endif
      enddo
      
      write(printstring,*) "Cannot find stream called ",tmpstr1(1:lgth)
      call printLine(WARNING,printstring)
      
      return
      end
