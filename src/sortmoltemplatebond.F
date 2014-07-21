      subroutine sortMolTemplateBond(atm1,res1,atm2,res2)
      implicit none
c#############################################################      
c
c     this routine alphabetically sorts at1 and at2, and
c     if they are the same, numerically sorts ires1,ires2
c
c     (C) Christopher Woods  18 February 2002
c
c##############################################################

      character*4 atm1,atm2,attemp,res1,res2,restemp
      
      if (atm1.lt.atm2) return
      if (atm1.eq.atm2) then
        if (res1.le.res2) return
      endif
      
      attemp = atm1
      atm1 = atm2
      atm2 = attemp
      
      restemp = res1
      res1 = res2
      res2 = restemp
      return
      end
