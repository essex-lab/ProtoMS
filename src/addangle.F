      subroutine addAngle(iprot,ires,tres,iang)
      include 'dimensions.inc'
      include 'templates.inc'
c############################################################
c
c     This routine adds a new angle within protein iprot,
c     residue ires, between atoms of type iang in template
c     tres...     
c
c     (C) Christopher Woods   23 January 2002
c
c############################################################

      integer iprot,tres,zat1,zat2,zat3
      integer ires,iang
      integer findAtom

c     get the zmatrix numbers of the three atoms in the angle...
c     if can't find the angle then return...
      zat1 = findAtom(iprot,ires,TemplateAngName(tres,iang,1))
      if (zat1.le.0) return
      zat2 = findAtom(iprot,ires,TemplateAngName(tres,iang,2))
      if (zat2.le.0) return
      zat3 = findAtom(iprot,ires,TemplateAngName(tres,iang,3))
      if (zat3.le.0) return

      call addProteinAngle(iprot,zat1,zat2,zat3,TemplateAngleType(tres,iang,1),
     .                TemplateAngleType(tres,iang,2),TemplateAngleDelta(tres,iang),
     .                TemplateAnglePar(tres,iang,1),TemplateAnglePar(tres,iang,2))

      return
      
      end
