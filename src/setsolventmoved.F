      subroutine setSolventMoved(isvn)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'movelist.inc'
      include 'printer.inc'
c###################################################
c
c     This sets the solvent 'ivn' as having moved.
c     This does *not* check to see if this solvent
c     has already been set as moved!
c
c     (C) Christopher Woods, October 2004
c
c###################################################

      integer isvn,in
      
c     do not need to do this if everything is already moved
      if (allMoved) return

      if (isvn.le.0 .or. isvn.gt.NSolvents) return
      
      in = nSvnMoved + 1
      if (in.gt.MAXSVNMOVED) then
        write(printstring,*) "Exceeded maximum number of moved solvents (",MAXSVNMOVED,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Increase MAXSVNMOVED and recompile")
        return
      endif
      
      svnMovedList(in) = isvn
      nSvnMoved = in
      return
      end
