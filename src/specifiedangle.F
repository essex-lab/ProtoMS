
      integer function specifiedAngle(ityp,ival,iat1,iat2,iat3)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'solutes.inc'
      include 'proteins.inc'
c########################################################################
c
c     This function returns whether there is a variable or 
c     additional angle specified between these three atoms, and
c     returns the index into the angle array for that angle
c
c     (C) Christopher Woods, January 2005
c
c########################################################################

       integer ityp,ival,iat1,iat2,iat3,nangs,i
       
       specifiedAngle = 0
       
       if (ityp.eq.PROTEIN) then
         nangs = NAngles(ival)
         do i=1,nangs
           if ( (Angle(ival,i,1).eq.iat1 .and. Angle(ival,i,2).eq.iat2 .and. Angle(ival,i,3).eq.iat3)
     .     .or. (Angle(ival,i,1).eq.iat3 .and. Angle(ival,i,2).eq.iat2 .and. Angle(ival,i,3).eq.iat1)) then
               
             specifiedAngle = i
             return
           endif
         enddo

         return
       else if (ityp.eq.SOLUTE) then
         nangs = NSolAngles(ival)
         do i=1,nangs
           if ( (SolAngle(ival,i,1).eq.iat1 .and. SolAngle(ival,i,2).eq.iat2
     .            .and. SolAngle(ival,i,3).eq.iat3)
     .     .or. (SolAngle(ival,i,1).eq.iat3 .and. SolAngle(ival,i,2).eq.iat2
     .            .and. SolAngle(ival,i,3).eq.iat1)) then
             
             specifiedAngle = i
             return
           endif
         enddo
       endif
       
       return
       
       end
