      logical function insolutesplitmoved(isol)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'quickpot.inc'
c########################################################
c
c     This returns true if solute isol is already
c     in the SplitMovedList
c
c     (C) Julien Michel , April 2005
c
c########################################################
      integer isol,i
            
      insolutesplitmoved = .false.

      if (nSolSplitMoved.le.0) then
        return
      else
        do i=1,nSolSplitMoved
          if (SplitSolMovedList(i).eq.isol)  then 
            insolutesplitmoved = .true.
            return
          endif
        enddo
c       only get here if the solute has not been moved
        return
      endif
      
      end
