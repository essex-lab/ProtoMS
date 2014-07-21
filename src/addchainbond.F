      subroutine addChainBond(iprot,ires,tres,ibnd)
      include 'dimensions.inc'
      include 'templates.inc'
c############################################################
c
c     This routine adds a new bond within protein iprot,
c     residue ires, between atoms of type ibnd in chain template
c     tres...     
c
c     (C) Christopher Woods
c
c############################################################

      integer iprot,tres,zat1,zat2
      integer ires,ibnd
      integer findAtom
      
c     get the zmatrix numbers of the two atoms in the bond...
c     if can't find the bond then return...
      zat1 = findAtom(iprot,ires,ChainBndName(tres,ibnd,1))
      if (zat1.le.0) return
      zat2 = findAtom(iprot,ires,ChainBndName(tres,ibnd,2))
      if (zat2.le.0) return

      call addProteinBond(iprot,zat1,zat2,ChainBondType(tres,ibnd,1),
     .             ChainBondType(tres,ibnd,2),ChainBondDelta(tres,ibnd),
     .             ChainBondPar(tres,ibnd,1),ChainBondPar(tres,ibnd,2))

      return
      
      end
