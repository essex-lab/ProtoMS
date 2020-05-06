      subroutine processidline()
      include 'dimensions.inc'
      include 'printer.inc'      
      include 'splitter.inc'
c#######################################
c     Read an id line and process it
c
c     Julien Michel February 2005
c
c######################################
      call lowercase(words(2))
      if (words(2) .eq. 'add') then
        call addid()
      else if (words(2) .eq. 'clear') then
        call clearid()
      else
        write(printstring,*) 'Could not interpret id line '
        call printLine(WARNING,printstring)
        call join(1,nwords,printstring)
        call printLine(WARNING,printstring)
      endif

      return
      end
