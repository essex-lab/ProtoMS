      subroutine genProteinAtom(iprot,zdih,zang,zbnd,iatz,idih,iang,ibnd)
      include 'dimensions.inc'
      include 'proteins.inc'
c##############################################################
c
c     This routine generates the coordinates of atom with 
c     zmatrix-id iatz, in protein iprot, based on the coordinate
c     of atoms with zmatrix id's zbnd, zang, zdih. 
c
c     iatz is bonded to zbnd, with bondlength ibnd A
c     iatz-zbnd-zang form an angle with size iang radians
c     iatz-zbnd-zang-zdih form a dihedral with size (counterclockwise) of 
c     idih radians.
c
c     (C) Christopher Woods   21 August 2002
c
c##############################################################

      integer iprot,zdih,zang,zbnd,iatz
      double precision ibnd,iang,idih
      double precision xbs,ybs,zbs,sinval
      double precision vecangdih(3),vecangbnd(3)
      double precision vx(3),vy(3),vz(3)
      integer i

c     zdih        iatm (generating coords of iatm!)
c        \        /
c       zang---zbnd

c     first create a set of x/y/z orthonormal vectors, with y perpendicular
c     vec(dih-bnd) and vec(ang-bnd), x perpendicular to vec(ang-bnd) and y, and z 
c     perpendicular to x and y. Do this via cross products...
      do i=1,3
        vecangdih(i) = ZCoords(iprot,zdih,i) - ZCoords(iprot,zbnd,i)
        vecangbnd(i) = ZCoords(iprot,zang,i) - ZCoords(iprot,zbnd,i)
      enddo

      call cross(vecangdih,vecangbnd,vy)
      call cross(vecangbnd,vy,vx)
      call cross(vx,vy,vz)
      
c     now we have the x/y/z vectors, we can generate the new coordinates
c     from this basis set...
c      thus x/y/z in the basis set is given by
c
c     xbs = bndlgth * sin(ang) * cos(dih)
c     ybs = bndlgth * sin(ang) * sin(dih)
c     zbs = -bndlgth * cos(ang)

      sinval = dsin(iang)
      xbs = ibnd * sinval * dcos(idih)
      ybs = ibnd * sinval * dsin(idih)
      zbs = -ibnd * dcos(iang)

c
c     Then we map the coordinates in this basis set to our cartesian coordinates
c     via...
c
c     x = xbs*vx(1) + ybs*vy(1) + zbs*vz(1)
c     y = xbs*vx(2) + ybs*vy(2) + zbs*vz(2)
c     z = xbs*vx(3) + ybs*vy(3) + zbs*vz(3)
c
c     These coordinates are based at the origin - they need to be based from 
c     the coordinates of the bond atom e.g.
c
c     x += xbnd     y+= ybnd   z+= zbnd
c
c     (we combine the last two steps together for speed)

      ZCoords(iprot,iatz,1) = vx(1)*xbs + vy(1)*ybs + vz(1)*zbs
     .                      + ZCoords(iprot,zbnd,1)
      ZCoords(iprot,iatz,2) = vx(2)*xbs + vy(2)*ybs + vz(2)*zbs
     .                      + ZCoords(iprot,zbnd,2)
      ZCoords(iprot,iatz,3) = vx(3)*xbs + vy(3)*ybs + vz(3)*zbs
     .                      + ZCoords(iprot,zbnd,3)
      
      return
      
      end
