      subroutine sortZDihedral(at1,at2,at3,at4)
      implicit none
      integer at1,at2,at3,at4,attemp
      
      if (at1.gt.at4) return
      if (at1.eq.at4) then
        if (at2.ge.at3) return
      endif
      
      attemp = at1
      at1 = at4
      at4 = attemp
      attemp = at2 
      at2 = at3
      at3 = attemp
      return
      end
