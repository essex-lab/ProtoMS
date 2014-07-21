      subroutine readSolTemplateCoordLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
c#############################################################
c
c     Unimplemented function - this will eventually allow
c     cartesian coordinates to be placed into the template
c     description
c
c#############################################################

      character*(*) line
      integer itmpl
      
      call printLine(WARNING,"readSolTemplateCoordLine() is not yet implemented!")
      write(printstring,*) "Arguments are line = ",line," and itmpl = ",itmpl

      return
      
      end
