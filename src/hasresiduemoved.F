      logical function hasResidueMoved(iprot,ires)
      implicit none
      include 'movelist.inc'
c########################################################
c
c     This returns true is residue ires in protein
c     iprot has been moved
c
c     (C) Christopher Woods, October 2004
c
c########################################################

      integer iprot,ires,i
            
      hasResidueMoved = .true.
      
      if (allMoved) return
      
      if (nResMoved.le.0) then
        hasResidueMoved = .false.
        return
      else
        do i=1,nResMoved
          if (resMovedList(i,1).eq.iprot .and. resMovedList(i,2).eq.ires) return
        enddo

c       only get here if the residue has not been moved
        hasResidueMoved = .false.
        return
      endif
      
      end
