      integer function findSA(type)
c########################################
c Return the index to the SA parameter 
c arrays for a given atom type. Return -1
c if the atom does not exist
c
c Julien Michel June 2004
c
c#########################################
      include "dimensions.inc"
      include "constants.inc"
      include "surface.inc"
      character*2 type
      integer i

      findSA = -1
      do 10 i=1,MAXSATYPES
        if (type .eq. SAType(i)) then
          findSA = i
          goto 11
c Faster look up by exiting if NULL character encountered
c this rely on the SAType array being filled sequentially
        else if (SAType(i) .eq. NULLCLJ) then 
          goto 11
        endif
 10   continue
 11   continue
      return 
      end
