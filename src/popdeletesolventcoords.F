      subroutine popDeleteSolventCoords
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'stacks.inc'
c#################################################
c
c     This deletes the top item on the coordinate
c     stack, and pops it
c
c     (C) Christopher Woods  9 June 2003
c
c#################################################

      if (svnCoordDepth.le.0) then
        call printLine(WARNING,"We have popped off the end of the solvent coordinate stack!")
        call printLine(WARNING,"This is a serious programming bug!")
        call closeProgram(-1,"ProtoMC is exiting as it has detected a serious programming bug!")
      endif
      
      svnCoordDepth = svnCoordDepth - 1
      return
      end
