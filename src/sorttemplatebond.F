      subroutine sortTemplateBond(at1,at2)
      implicit none
c###################################################
c
c     This subroutine alphabetically sorts at1 and at2
c
c     (C) Christopher Woods   22 January 2002
c
c###################################################

      character*4 at1,at2,at3
      
      if (at1.lt.at2) return
      at3 = at1
      at1 = at2
      at2 = at3
      return
      end
