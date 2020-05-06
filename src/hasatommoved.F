
      logical function hasAtomMoved(typ, id, atnum)
      include 'dimensions.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'movelist.inc'
      include 'enums.inc'
c################################################################################
c
c     This function returns whether or not an atom has moved. The atom
c     is identified via 'typ' (PROTEIN|SOLUTE|SOLVENT), 'id' (molecule number)
c     and 'atnum' (atom zmatrix number).
c
c     (C) Christopher Woods, November 2004
c
c################################################################################
      
      integer typ,id,atnum
      logical hasResidueMoved,hasSoluteMoved,hasSolventMoved

      if (allmoved) then
        hasAtomMoved = .true.
        return
      endif
      
      hasAtomMoved = .false.
      
      if (typ.eq.PROTEIN) then
        hasAtomMoved = hasResidueMoved(id,Zid2ResID(id,atnum,2))
      else if (typ.eq.SOLUTE) then
        hasAtomMoved = hasSoluteMoved(id)
      else if (typ.eq.SOLVENT) then
        hasAtomMoved = hasSolventMoved(id)
      endif
      
      return
      end
