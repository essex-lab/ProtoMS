      subroutine UpdateSplitMoveList(imov)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'quickpot.inc'
c##############################################
c     Save residues that have just moved 
c     if they haven't moved previously
c
c     (c) Julien Michel April 2005
c##############################################
      integer imov
c     Check what was the move done as this
      if (imov.eq.TRYPROTEIN) then
        call UpdateSplitMoveResidueList()
      else if (imov.eq.TRYSOLUTE) then
        call UpdateSplitMoveSoluteList()
      else if (imov.eq.TRYSOLVENT) then
        call UpdateSplitMoveSolventList()
      else
        write(printstring,*) "Problems sorting out the move type! ",imov
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is probably a programming bug!")
      endif

      return
      end
