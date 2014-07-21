
      integer function specifiedBond(ityp,ival,iat1,iat2)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'solutes.inc'
      include 'proteins.inc'
c########################################################################
c
c     This function returns whether there is a variable or 
c     additional bond specified between these two atoms, and
c     returns the index into the bond array for that bond
c
c     (C) Christopher Woods, January 2005
c
c########################################################################

       integer ityp,ival,iat1,iat2,nbnds,i
       
       specifiedBond = 0
       
       if (ityp.eq.PROTEIN) then
         nbnds = NBonds(ival)
         do i=1,nbnds
           if ( (Bond(ival,i,1).eq.iat1 .and. Bond(ival,i,2).eq.iat2)
     .         .or. (Bond(ival,i,1).eq.iat2 .and. Bond(ival,i,2).eq.iat1)) then
               
             specifiedBond = i
             return
           endif
         enddo

         return
       else if (ityp.eq.SOLUTE) then
         nbnds = NSolBonds(ival)
         do i=1,nbnds
           if ( (SolBond(ival,i,1).eq.iat1 .and. SolBond(ival,i,2).eq.iat2)
     .         .or. (SolBond(ival,i,1).eq.iat2 .and. SolBond(ival,i,2).eq.iat1)) then
             
             specifiedBond = i
             return
           endif
         enddo
       endif
       
       return
       
       end
