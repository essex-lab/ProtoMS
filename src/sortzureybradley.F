      subroutine sortZUreyBradley(at1,at2,at3)
      implicit none
      integer at1,at2,at3,attemp
      
      if (at1.ge.at3) return
      attemp = at1
      at1 = at3
      at3 = attemp
      return
      end
