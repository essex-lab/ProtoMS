
      character*6 function positionString(ipos)
      implicit none
      include 'enums.inc'
c####################################################
c
c     This function returns a string describing
c     the position ipos in the protein chain
c     (e.g. 'first', 'middle' etc.)
c
c     (C) Christopher Woods, November 2004
c
c####################################################

      integer ipos
      
      if (ipos.eq.FIRST) then
        positionString = "first"
      else if (ipos.eq.MIDDLE) then
        positionString = "middle"
      else if (ipos.eq.LAST) then
        positionString = "last"
      else if (ipos.eq.SINGLE) then
        positionString = "single"
      else if (ipos.eq.ALLPOS) then
        positionString = "all"
      else
        positionString = "none"
      endif
      
      return
      end
