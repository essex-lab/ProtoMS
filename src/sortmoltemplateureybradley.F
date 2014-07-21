      subroutine sortMolTemplateUreyBradley(atm1,res1,atm2,res2,atm3,res3)
      implicit none
c#############################################################      
c
c     this routine alphabetically sorts at1 and at2, and
c     if they are the same, numerically sorts ires1,ires2
c
c     (C) Christopher Woods
c
c##############################################################

      character*4 res1,res2,res3,restemp
      character*4 atm1,atm2,atm3,atmtemp
      
      if (atm1.lt.atm3) return
      if (atm1.eq.atm3) then
        if (res1.le.res3) return
      endif
      
      atmtemp = atm1
      atm1 = atm3
      atm3 = atmtemp
      
      restemp = res1
      res1 = res3
      res3 = restemp
      return
      end
