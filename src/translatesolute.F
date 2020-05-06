      subroutine translateSolute(isol,idx)
      include 'dimensions.inc'
      include 'solutes.inc'
      integer isol,i,j
      double precision idx(3)

      do j=1,3
        do i=1,3
          SolDummy(isol,i,j) = SolDummy(isol,i,j) + idx(j)
        enddo
        SolGeomCentre(isol,j) = SolGeomCentre(isol,j) + idx(j)
      enddo 

      return
      end
