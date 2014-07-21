      subroutine sortTemplateDihedral(at1,at2,at3,at4)
      implicit none
c###################################################
c
c     This subroutine alphabetically sorts the dihedral
c     angle at1,at2,at3,at4
c
c     If at1<at4 leave alone
c     if at1>at4 then reverse order
c     if (at1=at4) then do the same test for at2 and at3
c
c     (C) Christopher Woods   22 January 2002
c
c###################################################

      character*4 at1,at2,at3,at4,at5
      
      if (at1.lt.at4) return
      if (at1.eq.at4) then
        if (at2.le.at3) return
      endif
      
c     we have to swap the angle...
c     swap the ends first...
      at5 = at1
      at1 = at4
      at4 = at5
c     now swap the middle
      at5 = at2
      at2 = at3
      at3 = at5
      return

      end
