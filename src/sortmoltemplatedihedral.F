      subroutine sortMolTemplateDihedral(atm1,res1,atm2,res2,atm3,res3,atm4,res4)
      implicit none
c################################################################
c
c     Dihedral equivelent of sortMolTemplateBond
c
c     (C) Christopher Woods  18 February 2002
c
c################################################################

      character*4 atm1,atm2,atm3,atm4,atmtemp
      character*4 res1,res2,res3,res4,restemp
      
      if (atm1.lt.atm4) return
      if (atm1.eq.atm4) then
        if (res1.lt.res4) return
        if (res1.eq.res4) then
          if (atm2.lt.atm3) return
            if (atm2.eq.atm3) then
              if (res2.le.res3) return
          endif
        endif
      endif
      
      atmtemp = atm1
      atm1 = atm4
      atm4 = atmtemp
      atmtemp = atm2 
      atm2 = atm3
      atm3 = atmtemp
      
      restemp = res1
      res1 = res4
      res4 = restemp
      restemp = res2
      res2 = res3
      res3 = restemp
      return
      end
