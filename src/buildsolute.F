      subroutine buildSolute(isol)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'solutes.inc'
c####################################################
c
c     This subroutine rebuilds the solute isol
c     from the zmatrix internal coordinates and
c     the position of the calculated dummy atoms...
c
c     (C) Christopher Woods   28 February 2002
c
c####################################################

      integer isol
      integer i,nats

      if (isol.le.0 .or. isol.gt.MAXSOLUTES) return
      
      nats = NSolAtoms(isol)
      if (nats.le.0) return
      
      do i=1,nats
        call genSoluteAtom(isol,ZSolMat(isol,i,DIHEDRALS),ZSolMat(isol,i,ANGLES),
     .                          ZSolMat(isol,i,BONDS),i,
     .                          ZSolVal(isol,i,DIHEDRALS),ZSolVal(isol,i,ANGLES),
     .                          ZSolVal(isol,i,BONDS))
      enddo

      call updateSoluteInternals(isol)

c     get the new centre and radius of the solute
      call getSoluteCentreRadius(isol)
      
      return
      end
