      subroutine addUreyBradley(iprot,ires,tres,iuby)
      include 'dimensions.inc'
      include 'templates.inc'
c############################################################
c
c     This routine adds a new urey-bradley within protein iprot,
c     residue ires, between atoms of type iuby in template
c     tres...     
c
c     (C) Christopher Woods
c
c############################################################

      integer iprot,tres,zat1,zat2,zat3
      integer ires,iuby
      integer findAtom
      
c     get the zmatrix numbers of the two atoms in the bond...
c     if can't find the bond then return...
      zat1 = findAtom(iprot,ires,TemplateUByName(tres,iuby,1))
      if (zat1.le.0) return
      zat2 = findAtom(iprot,ires,TemplateUByName(tres,iuby,2))
      if (zat2.le.0) return
      zat3 = findAtom(iprot,ires,TemplateUByName(tres,iuby,3))
      if (zat3.le.0) return

      call addProteinUreyBradley(iprot,zat1,zat2,zat3,TemplateUreyBradleyType(tres,iuby),
     .          TemplateUreyBradleyPar(tres,iuby,1),TemplateUreyBradleyPar(tres,iuby,2))

      return
      
      end
