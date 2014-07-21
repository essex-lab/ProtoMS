      logical function hasSolventMoved(isvn)
      implicit none
      include 'movelist.inc'
c########################################################
c
c     This returns true is solvent isvn has been moved
c
c     (C) Christopher Woods, October 2004
c
c########################################################

      integer isvn,i
            
      hasSolventMoved = .true.
      
      if (allMoved) return
      
      if (nSvnMoved.le.0) then
        hasSolventMoved = .false.
        return
      else
        do i=1,nSvnMoved
          if (svnMovedList(i).eq.isvn) return
        enddo

c       only get here if the solvent has not been moved
        hasSolventMoved = .false.
        return
      endif
      
      end
