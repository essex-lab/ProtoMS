
      subroutine setProteinGeometry(iprot,ires,iatmnam,ityp,ival)
      include 'dimensions.inc'
      include 'proteins.inc'
      include 'printer.inc'
      include 'enums.inc'
c######################################################################
c
c     Set a bond/angle/dihedral within protein 'iprot' to the value
c     'ival'. ires is the residue number of the atom, iatmnam,
c     that will have its zmatrix geometry (type ityp) set the value
c     of ival. If ityp is BONDS, then ival is in angstroms, while
c     if ityp is ANGLES or DIHEDRALS then ival is in degrees.
c
c     (C) Christopher Woods, January 2005
c
c######################################################################

      integer iprot,ityp,iatm,ires
      integer findAtom,findResidue
      character*4 iatmnam
      double precision ival,deg2Rad
      
c     get the index into the residue array for residue numbered ires
      ires = findResidue(iprot,ires)
      if (ires.le.0) then
        write(printstring,*) "Could not find protein residue ",ires
        call printLine(WARNING,printstring)
        return
      endif
      
c     get the zmatrix atom number for this atom...
      iatm = findAtom(iprot,ires,iatmnam)
      
      if (iatm.eq.0) then
        write(printstring,*) "Could not find protein atom ",iatmnam," - ",ires
        call printLine(WARNING,printstring)
        return
      endif

c     now get the index of this atom in the zmatrix arrays
      iatm = Zid2ResID(iprot,iatm,2)
      if (iatm.le.0) then
        call printLine(WARNING,"Can only set the geometry of sidechain atoms!")
        return
      endif
      
      if (ityp.eq.BONDS) then
        ZValSC(iprot,ires,iatm,BONDS) = ival        
      else if (ityp.eq.ANGLES .or. ityp.eq.DIHEDRALS) then
        ZValSC(iprot,ires,iatm,ityp) = deg2Rad(ival)
      endif

c     now rebuild the coordinates of the solute      
      call buildResidue(iprot,ires)           

      return
      
      end
