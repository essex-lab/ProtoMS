      subroutine getVectorOnSphere(vec)
c     This is taken from Frenkel Smitt p578
c
c     Julien Michel October 2004
      double precision vec(3)
      double precision ransq,ran1,ran2,ranh
      double precision ranRange
      ransq = 2.0D+0
      ran1 = 0.0D+0
      ran2 = 0.0D+0
      do while (ransq .ge. 1.0D+0)
        ran1 = 1 - 2.0D+0 * ranRange(0.0D+0,1.0D+0)
        ran2 = 1 - 2.0D+0 * ranRange(0.0D+0,1.0D+0)
        ransq = ran1 ** 2 + ran2 ** 2
      enddo
      ranh = 2 * sqrt( 1.0D+0 - ransq )
      vec(1) = ran1 * ranh
      vec(2) = ran2 * ranh
      vec(3) = (1 - 2.0D+0*ransq)

      return
      end
