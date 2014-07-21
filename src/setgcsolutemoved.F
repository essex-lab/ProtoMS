      subroutine setGCSoluteMoved(igcs)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'movelist.inc'
      include 'printer.inc'
c###################################################
c
c     This sets the GC-solute 'igcs' as having moved.
c     This does *not* check to see if this solvent
c     has already been set as moved!
c
c     (C) Based on routine of Christopher Woods, October 2004
c      MSB August 2012
c###################################################

      integer igcs,in
      
c     do not need to do this if everything is already moved
      if (allMoved) return

      if (igcs.le.0 .or. igcs.gt.NGCSolutes) return
      
      in = nGCSMoved + 1
      if (in.gt.MAXSVNMOVED) then
        write(printstring,*) "Exceeded maximum number of moved GC-solutes (",MAXSVNMOVED,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Increase MAXSVNMOVED and recompile")
        return
      endif
      
      gcsMovedList(in) = igcs
      nGCSMoved = in
      return
      end
