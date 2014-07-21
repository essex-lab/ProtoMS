      subroutine addBond(iprot,ires,tres,ibnd)
      include 'dimensions.inc'
      include 'templates.inc'
c############################################################
c
c     This routine adds a new bond within protein iprot,
c     residue ires, between atoms of type ibnd in template
c     tres...     
c
c     (C) Christopher Woods   23 January 2002
c
c############################################################

      integer iprot,tres,zat1,zat2
      integer ires,ibnd
      integer findAtom
      
c     get the zmatrix numbers of the two atoms in the bond...
c     if can't find the bond then return...
      zat1 = findAtom(iprot,ires,TemplateBndName(tres,ibnd,1))
      if (zat1.le.0) return
      zat2 = findAtom(iprot,ires,TemplateBndName(tres,ibnd,2))
      if (zat2.le.0) return

      call addProteinBond(iprot,zat1,zat2,TemplateBondType(tres,ibnd,1),
     .             TemplateBondType(tres,ibnd,2),TemplateBondDelta(tres,ibnd),
     .             TemplateBondPar(tres,ibnd,1),TemplateBondPar(tres,ibnd,2))

      return
      
      end
