      integer function findGB(type)
c########################################
c Return the index to the GB parameter 
c arrays for a given atom type. Return -1
c if the atom does not exist
c
c Julien Michel November 2004
c
c#########################################
      include "dimensions.inc"
      include "constants.inc"
      include "gb.inc"
      character*2 type
      integer i

      findGB = -1
      do 10 i=1,MAXGBTYPES
        if (type .eq. GBType(i)) then
          findGB = i
          goto 11
c Faster look up by exiting if NULL character encountered
c this rely on the SAType array being filled sequentially
        else if (GBType(i) .eq. NULLCLJ) then 
          goto 11
        endif
 10   continue
 11   continue
      return 
      end
