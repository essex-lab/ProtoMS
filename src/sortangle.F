      subroutine sortAngle(at1,at2,at3)
      implicit none
c###################################################
c
c     This subroutine alphabetically sorts at1 and at3,
c     leaving at2 in the middle
c
c     (C) Christopher Woods   22 January 2002
c
c###################################################

      character*2 at1,at2,at3,at4
      
      if (at1.lt.at3) return
      at4 = at1
      at1 = at3
      at3 = at4
      return
      end
