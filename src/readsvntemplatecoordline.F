      subroutine readSvnTemplateCoordLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
c############################################################
c
c     This will eventually allow the reading in of 
c     cartesian coordinates for the template solvent
c
c############################################################

      character*(*) line
      integer itmpl
      
      call printLine(WARNING,
     .      "readSvnTemplateCoordLine() is not yet implemented - would you like to do it?")
     
      return
      end
