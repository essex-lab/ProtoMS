      subroutine getSoluteZmatrix(isol)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'solutes.inc'
c##########################################################
c
c     This routine gets the zmatrix values (bonds/angles/dihedrals)
c     for the current solute, from its current pdb coordinates.
c
c     (C) Christopher Woods   28 February 2002
c
c##########################################################

      integer isol,i,nats
      double precision calcBond,calcAngle,calcDihedral

      if (isol.le.0 .or. isol.gt.MAXSOLUTES) return
      
      nats = NSolAtoms(isol)
      if (nats.le.0) return

      do i=1,nats
        ZSolVal(isol,i,1) = calcBond(SOLUTE,isol,i,ZSolMat(isol,i,1))
        ZSolVal(isol,i,2) = calcAngle(SOLUTE,isol,i,ZSolMat(isol,i,1),ZSolMat(isol,i,2))
        ZSolVal(isol,i,3) = calcDihedral(SOLUTE,isol,i,ZSolMat(isol,i,1),
     .                ZSolMat(isol,i,2),ZSolMat(isol,i,3))
      enddo

      return
      end
