      subroutine setLength(ilgth,vec)
      include 'dimensions.inc'
      include 'constants.inc'
      double precision ilgth,ilgth2,vec(3)
      double precision veclength
      integer i
c##########################################

      ilgth2 = veclength(vec)
      
      if (ilgth2.le.ZERO) return
      
      ilgth2 = ilgth / ilgth2
      
      do i=1,3
        vec(i) = ilgth2 * vec(i)
      enddo

      return
      end
