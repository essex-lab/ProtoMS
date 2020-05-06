      double precision function veclength2(vec)
      include 'dimensions.inc'
      include 'constants.inc'
      double precision vec(3)
      double precision length
      integer i
c##########################################

      length = ZERO
      do i=1,3
        length = length + (vec(i)**2)
      enddo

      veclength2 = length
      
      end
