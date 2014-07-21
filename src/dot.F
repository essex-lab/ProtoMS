      double precision function dot(a,b)
      implicit none
c###############################################################
c
c     Calculates the dot product of the two vectors, a and b
c
c     (C) Christopher Woods       21 August 2002
c
c###############################################################

      double precision a(3),b(3)
      double precision temp
      integer i

      temp = 0.D+00
      do i=1,3
        temp = temp + (a(i)*b(i))
      enddo

      dot = temp
      
      end
