      subroutine setSoluteMoved(isol)
      include 'dimensions.inc'
      include 'solutes.inc'
      include 'movelist.inc'
      include 'printer.inc'
c###################################################
c
c     This sets the solute 'isol' as having moved.
c     This does *not* check to see if this solute
c     has already been set as moved!
c
c     (C) Christopher Woods, October 2004
c
c###################################################

      integer isol,in
      
c     do not need to do this if everything is already moved
      if (allMoved) return

      if (isol.le.0 .or. isol.gt.NSolutes) return
      
      in = nSolMoved + 1
      if (in.gt.MAXSOLMOVED) then
        write(printstring,*) "Exceeded maximum number of moved solutes (",MAXSOLMOVED,")"
        call printLine(FATAL,printstring)
        write(printstring,*)"Increase MAXSOLMOVED and recompile"
        call closeProgram(-1,printstring)        
        return
      endif
      
      solMovedList(in) = isol
      nSolMoved = in
      return
      end
