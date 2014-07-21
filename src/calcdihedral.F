      double precision function calcDihedral(ityp,ival,izat1,izat2,izat3,izat4)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'printer.inc'
c##############################################################################
c
c     This function returns the size of the dihedral angle (counterclockwise, radians)
c    
c     (C) Christopher Woods
c
c##############################################################################

      integer iprot,izat1,izat2,izat3,izat4,i,ityp,ival,isol
      double precision at1(3),at2(3),at3(3),at4(3)
      double precision calcGeneralDihedral
      
c     first, copy the coords of the atoms into at1,at2,at3,at4      
      if (ityp.eq.PROTEIN) then
        iprot = ival
        do i=1,3
          at1(i) = ZCoords(iprot,izat1,i)
          at2(i) = ZCoords(iprot,izat2,i)
          at3(i) = ZCoords(iprot,izat3,i)
          at4(i) = ZCoords(iprot,izat4,i)
        enddo
      else if (ityp.eq.SOLUTE) then
        isol = ival
        if (izat1.lt.0) then
          do i=1,3
            at1(i) = SolDummy(isol,-izat1,i)
          enddo
        else
          do i=1,3
            at1(i) = ZSolCoords(isol,izat1,i)
          enddo
        endif
        if (izat2.lt.0) then
          do i=1,3
            at2(i) = SolDummy(isol,-izat2,i)
          enddo
        else
          do i=1,3
            at2(i) = ZSolCoords(isol,izat2,i)
          enddo
        endif
        if (izat3.lt.0) then
          do i=1,3
            at3(i) = SolDummy(isol,-izat3,i)
          enddo
        else
          do i=1,3
            at3(i) = ZSolCoords(isol,izat3,i)
          enddo
        endif
        if (izat4.lt.0) then
          do i=1,3
            at4(i) = SolDummy(isol,-izat4,i)
          enddo
        else
          do i=1,3
            at4(i) = ZSolCoords(isol,izat4,i)
          enddo
        endif
      else
        write(printstring,*) "Unknown system requested in calcDihedral() - ",ityp
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is probably a programming bug!")
        calcDihedral = ZERO
        return
      endif

      calcDihedral = calcGeneralDihedral(at1,at2,at3,at4)

c     JM Feb 2005 . The code below has been moved to calcgeneraldihedral()
c$$$c     
c$$$c     at1        at4
c$$$c       \        /
c$$$c       at2----at3  
c$$$c
c$$$c     Dihedral angle is the plane of intersection between the planes
c$$$c     formed by at1,at2,at3 and at2,at3,at4.
c$$$c
c$$$c     This is equivelent to the angle of intersection of the normals to those
c$$$c     planes - thus we must calculate those normals!
c$$$c
c$$$c     normal_1 = cross(vec(at2->at1), vec(at2->at3))
c$$$c     normal_2 = cross(vec(at3->at2), vec(at3->at4))
c$$$      
c$$$      do i=1,3
c$$$        vec21(i) = at1(i) - at2(i)
c$$$        vec23(i) = at3(i) - at2(i)
c$$$        vec32(i) = -vec23(i)
c$$$        vec34(i) = at4(i) - at3(i)
c$$$      enddo
c$$$
c$$$      call cross(vec21,vec23,norm1)
c$$$      call cross(vec32,vec34,norm2)
c$$$
c$$$c     now the phi is the angle of intersection of the normals, e.g.
c$$$c
c$$$c     cos(phi) = dot(normal_1,normal_2) / (length(normal_1)*length(normal_2))
c$$$c       (and length of each normal is 1.0)
c$$$      
c$$$      phi = dot(norm1,norm2)
c$$$
c$$$c     need to take the inv_cos of cos(phi) - thus must make sure that it
c$$$c     lies between -1 and 1                        
c$$$      if (phi.gt.ONE) then
c$$$        phi = ONE
c$$$      else if (phi.lt.-ONE) then
c$$$        phi = -ONE
c$$$      endif
c$$$        
c$$$      phi = dacos(phi)
c$$$
c$$$c     this only tells us if the angle lies between 0 and pi (0 and 180) 
c$$$c     To get the full size of the angle (whether > pi, >180) we need to calculate
c$$$c     the angle vs. the norm of norm_1 and vec23, and see if this is less than pi/2 (90)
c$$$c     If it is, then the angle is >pi, <2pi (>180, <360). 
c$$$      call cross(norm1,vec23,norm3)
c$$$      ang = dot(norm2,norm3)
c$$$      if (ang.gt.ONE) then
c$$$        ang = ONE
c$$$      else if (ang.lt.-ONE) then
c$$$        ang = -ONE
c$$$      endif
c$$$      ang = dacos(ang)
c$$$      
c$$$      if (ang.lt.PIOTWO) phi = TWOPI-phi
c$$$                                               
c$$$      calcDihedral = phi
      
      return
      
      end
