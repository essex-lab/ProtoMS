      logical function hasSoluteMoved(isol)
      implicit none
      include 'movelist.inc'
c########################################################
c
c     This returns true is solute isol has been moved
c
c     (C) Christopher Woods, October 2004
c
c########################################################

      integer isol,i
            
      hasSoluteMoved = .true.
      
      if (allMoved) return
      
      if (nSolMoved.le.0) then
        hasSoluteMoved = .false.
        return
      else
        do i=1,nSolMoved
          if (solMovedList(i).eq.isol) return
        enddo

c       only get here if the solute has not been moved
        hasSoluteMoved = .false.
        return
      endif
      
      end
