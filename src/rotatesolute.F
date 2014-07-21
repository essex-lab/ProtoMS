      subroutine rotateSolute(isol,ranvec,ranval)
      include 'dimensions.inc'
      include 'solutes.inc'
      integer isol,i,j
      double precision ranvec(3),ranval,ivec(3)
      double precision vecLength,ilgth

      do j=1,3
        do i=1,3
          ivec(i) = SolDummy(isol,j,i) - SolGeomCentre(isol,i)
        enddo
        ilgth = vecLength(ivec)
        call rotateVector(ivec,ranvec,ranval)
c       restore the length of the vector - prevents precision
c       errors from propagating
        call setLength(ilgth,ivec)
        do i=1,3
          SolDummy(isol,j,i) = SolGeomCentre(isol,i) + ivec(i)
        enddo
      enddo
            
      return
      end
