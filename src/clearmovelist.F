      subroutine clearMoveList()
      include 'dimensions.inc'
      include 'movelist.inc'
c################################################
c
c     This routine clears the move lists
c
c     (C) Christopher Woods
c
c################################################

      nResMoved = 0
      nSolMoved = 0
      nSvnMoved = 0
      nGCSMoved = 0
      allmoved = .false.
      
      return
      end
