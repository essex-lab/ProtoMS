      subroutine sortMolTemplateAngle(atm1,res1,atm2,res2,atm3,res3)
      implicit none
c##############################################################
c
c     Angle equivelent of sortMolTemplateBond
c
c     (C) Christopher Woods   18 February 2002
c
c##############################################################

      character*4 atm1,atm2,atm3,attemp
      character*4 res1,res2,res3,restemp
        
      if (atm1.lt.atm3) return
      if (atm1.eq.atm3) then
        if (res1.le.res3) return
      endif
      
      attemp = atm1
      atm1 = atm3
      atm3 = attemp
      
      restemp = res1
      res1 = res3
      res3 = restemp
      return
      end
