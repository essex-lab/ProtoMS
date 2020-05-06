      double precision function calcGeneralDihedral(at1,at2,at3,at4)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'constants.inc'
      include 'printer.inc'
c##############################################################################
c
c     This function returns the size of the dihedral angle 
c     (counterclockwise, radians)
c     Unlike calcdihedral, the atoms forming the dihedral do not have to 
c     be in the same molecule.
c    
c     Julien Michel February 2005, taken from calcdihedral
c
c##############################################################################
      integer i
      double precision at1(3),at2(3),at3(3),at4(3)
      double precision vec21(3),vec23(3),vec32(3),vec34(3)
      double precision norm1(3),norm2(3),norm3(3),phi,ang
      double precision dot

c
c     at1        at4
c       \        /
c       at2----at3  
c
c     Dihedral angle is the plane of intersection between the planes
c     formed by at1,at2,at3 and at2,at3,at4.
c
c     This is equivelent to the angle of intersection of the normals to those
c     planes - thus we must calculate those normals!
c
c     normal_1 = cross(vec(at2->at1), vec(at2->at3))
c     normal_2 = cross(vec(at3->at2), vec(at3->at4))
      
      do i=1,3
        vec21(i) = at1(i) - at2(i)
        vec23(i) = at3(i) - at2(i)
        vec32(i) = -vec23(i)
        vec34(i) = at4(i) - at3(i)
      enddo

      call cross(vec21,vec23,norm1)
      call cross(vec32,vec34,norm2)

c     now the phi is the angle of intersection of the normals, e.g.
c
c     cos(phi) = dot(normal_1,normal_2) / (length(normal_1)*length(normal_2))
c       (and length of each normal is 1.0)
      
      phi = dot(norm1,norm2)

c     need to take the inv_cos of cos(phi) - thus must make sure that it
c     lies between -1 and 1                        
      if (phi.gt.ONE) then
        phi = ONE
      else if (phi.lt.-ONE) then
        phi = -ONE
      endif
        
      phi = dacos(phi)

c     this only tells us if the angle lies between 0 and pi (0 and 180) 
c     To get the full size of the angle (whether > pi, >180) we need to calculate
c     the angle vs. the norm of norm_1 and vec23, and see if this is less than pi/2 (90)
c     If it is, then the angle is >pi, <2pi (>180, <360). 
      call cross(norm1,vec23,norm3)
      ang = dot(norm2,norm3)
      if (ang.gt.ONE) then
        ang = ONE
      else if (ang.lt.-ONE) then
        ang = -ONE
      endif
      ang = dacos(ang)
      
      if (ang.lt.PIOTWO) phi = TWOPI-phi
                                               
      calcGeneralDihedral = phi
      
      return
      
      end
