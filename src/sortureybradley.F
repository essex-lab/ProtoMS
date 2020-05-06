      subroutine sortUreyBradley(at1,at2,at3)
      implicit none
c###################################################
c
c     This subroutine alphabetically sorts at1 and at2
c
c     (C) Christopher Woods
c
c###################################################

      character*4 at1,at2,at3,at4
      
      if (at1.lt.at3) return
      at4 = at1
      at1 = at3
      at3 = at4
      return
      end
