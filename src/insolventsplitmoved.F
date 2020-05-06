      logical function insolventsplitmoved(isvn)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'quickpot.inc'
c########################################################
c
c     This returns true if solvent isvn is already
c     in the SplitMovedList
c
c     (C) Julien Michel , April 2005
c
c########################################################
      integer isvn,i
            
      insolventsplitmoved = .false.

      if (nSvnSplitMoved.le.0) then
        return
      else
        do i=1,nSvnSplitMoved
          if (SplitSvnMovedList(i).eq.isvn)  then 
            insolventsplitmoved = .true.
            return
          endif
        enddo
c       only get here if the solvent has not been moved
        return
      endif
      
      end
