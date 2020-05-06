      subroutine sortZBond(at1,at2)
      implicit none
      integer at1,at2,attemp
      
      if (at1.ge.at2) return
      attemp = at1
      at1 = at2
      at2 = attemp
      return
      end
