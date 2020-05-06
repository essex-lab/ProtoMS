      subroutine prettyLine(stream,line)
      include 'dimensions.inc'
      include 'printer.inc'
c##################################################################
c
c     This routine does the same as printLine, except it 
c     draws a pretty box around the line to highlight it!
c
c       (C) Christopher Woods, October 2004
c
c##################################################################

      integer stream
      character*(*) line
      character*300 prettystring
      integer lineLength,FILE
      integer lgth,i

      if (.not.prettyprint) then
c       just do normal printing if pretty printing is disabled        
        call printLine(stream,line)
        return
      endif
      
c     get the file associated with this stream
      FILE = streamhandle(stream)
      
c     do nothing if there is no file!
      if (FILE.lt.0) return
            
c     now get the length of text in printstring
      lgth = lineLength(line)
      
      if (lgth.le.1) then
        return
      else if (lgth.gt.300) then
        lgth = 300
      endif
      
      write(prettystring,11) ('#',i=1,lgth)
11    format(1000a1)

      call printLine(stream,prettystring)
      call printLine(stream,line)
      call printLine(stream,prettystring)
  
      end
