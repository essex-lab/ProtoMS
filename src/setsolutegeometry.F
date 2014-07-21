
      subroutine setSoluteGeometry(isol,iresnam,iatmnam,ityp,ival)
      include 'dimensions.inc'
      include 'solutes.inc'
      include 'printer.inc'
      include 'enums.inc'
c######################################################################
c
c     Set a bond/angle/dihedral within solute 'isol' to the value
c     'ival'. iresnam is the residue name of the atom, iatmnam,
c     that will have its zmatrix geometry (type ityp) set the value
c     of ival. If ityp is BONDS, then ival is in angstroms, while
c     if ityp is ANGLES or DIHEDRALS then ival is in degrees.
c
c     (C) Christopher Woods, January 2005
c
c######################################################################

      integer isol,ityp,iatm
      integer findSoluteAtom
      character*4 iresnam,iatmnam
      double precision ival,deg2Rad
      
c     get the zmatrix atom number for this atom...
      iatm = findSoluteAtom(isol,iatmnam,iresnam)
      if (iatm.eq.0) then
        write(printstring,*) "Could not find solute atom ",iatmnam," - ",iresnam
        call printLine(WARNING,printstring)
        return
      else if (iatm.lt.0) then
        write(printstring,*) "Cannot change geometry of solute dummy atoms!"
        call printLine(WARNING,printstring)
        return
      endif

      if (ityp.eq.BONDS) then
        ZSolVal(isol,iatm,BONDS) = ival        
      else if (ityp.eq.ANGLES .or. ityp.eq.DIHEDRALS) then
        write(printstring,*) "Setting geometry to ",ival," degrees, ",deg2Rad(ival)," radians"
        call printLine(INFO,printstring)
        ZSolVal(isol,iatm,ityp) = deg2Rad(ival)
      endif

c     now rebuild the coordinates of the solute      
      call buildSolute(isol)            

      return
      
      end
