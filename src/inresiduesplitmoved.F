      logical function inresiduesplitmoved(iprot,ires)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'quickpot.inc'
c########################################################
c
c     This returns true if residue ires of iprot is already
c     in the SplitMovedList
c
c     (C) Julien Michel , April 2005
c
c########################################################
      integer iprot,i,ires
      
      inresiduesplitmoved = .false.

      if (nResSplitMoved.le.0) then
        return
      else
        do i=1,nResSplitMoved
          if (SplitResMovedList(i,1).eq.iprot .and. 
     .      SplitResMovedList(i,2) .eq. ires)  then 
            inresiduesplitmoved = .true.
            return
          endif
        enddo
c       only get here if the residue has not been moved
        return
      endif
      
      end
