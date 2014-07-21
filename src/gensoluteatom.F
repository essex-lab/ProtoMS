      subroutine genSoluteAtom(isol,zdih,zang,zbnd,iatz,idih,iang,ibnd)
      include 'dimensions.inc'
      include 'solutes.inc'
c#####################################################################
c
c     This routine generates the coordinates of atom with zmatrix
c     id iatz, in solute isol, based on coordinates of atoms
c     zbnd,zang,zdih, and bond/angle/dihedral ibnd,iang,idih
c
c     (C) Christopher Woods 28 February 2002
c
c#####################################################################

      integer isol,zdih,zang,zbnd,iatz
      double precision ibnd,iang,idih
      double precision xbs,ybs,zbs,sinval
      double precision vecangdih(3),vecangbnd(3)
      double precision vx(3),vy(3),vz(3)
      double precision ia(3),ib(3),ic(3)
      integer i

c     zdih        iatm (generating coords of iatm!)
c        \        /
c       zang---zbnd

c     first create a set of x/y/z orthonormal vectors, with y perpendicular
c     vec(dih-bnd) and vec(ang-bnd), x perpendicular to vec(ang-bnd) and y, and z 
c     perpendicular to x and y. Do this via cross products...
      if (zbnd.lt.0) then
        do i=1,3
          ia(i) = SolDummy(isol,-zbnd,i)
        enddo
      else
        do i=1,3
          ia(i) = ZSolCoords(isol,zbnd,i)
        enddo
      endif
      
      if (zang.lt.0) then
        do i=1,3
          ib(i) = SolDummy(isol,-zang,i)
        enddo
      else
        do i=1,3
          ib(i) = ZSolCoords(isol,zang,i)
        enddo
      endif
      
      if (zdih.lt.0) then
        do i=1,3
          ic(i) = SolDummy(isol,-zdih,i)
        enddo
      else
        do i=1,3
          ic(i) = ZSolCoords(isol,zdih,i)
        enddo
      endif
      
      do i=1,3
        vecangdih(i) = ic(i) - ia(i)
        vecangbnd(i) = ib(i) - ia(i)
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

      ZSolCoords(isol,iatz,1) = vx(1)*xbs + vy(1)*ybs + vz(1)*zbs
     .                      + ia(1)
      ZSolCoords(isol,iatz,2) = vx(2)*xbs + vy(2)*ybs + vz(2)*zbs
     .                      + ia(2)
      ZSolCoords(isol,iatz,3) = vx(3)*xbs + vy(3)*ybs + vz(3)*zbs
     .                      + ia(3)
      
      return
      
      end
