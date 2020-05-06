      subroutine addChainAngle(iprot,ires,tres,iang)
      include 'dimensions.inc'
      include 'templates.inc'
c############################################################
c
c     This routine adds a new angle within protein iprot,
c     residue ires, between atoms of type iang in chain template
c     tres...     
c
c     (C) Christopher Woods
c
c############################################################

      integer iprot,tres,zat1,zat2,zat3
      integer ires,iang
      integer findAtom

c     get the zmatrix numbers of the three atoms in the angle...
c     if can't find the angle then return...
      zat1 = findAtom(iprot,ires,ChainAngName(tres,iang,1))
      if (zat1.le.0) return
      zat2 = findAtom(iprot,ires,ChainAngName(tres,iang,2))
      if (zat2.le.0) return
      zat3 = findAtom(iprot,ires,ChainAngName(tres,iang,3))
      if (zat3.le.0) return

      call addProteinAngle(iprot,zat1,zat2,zat3,ChainAngleType(tres,iang,1),
     .                ChainAngleType(tres,iang,2),ChainAngleDelta(tres,iang),
     .                ChainAnglePar(tres,iang,1),ChainAnglePar(tres,iang,2))

      return
      
      end
