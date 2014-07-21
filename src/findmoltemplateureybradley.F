      integer function findMolTemplateUreyBradley(atm1,res1,atm2,res2,atm3,res3,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c########################################################
c
c     This finds the Urey-Bradley term between at1(res1)-at2(res2) in
c     molecule template itmpl, and returns the index number,
c     or 0 if it does not exist
c
c     (C) Christopher Woods
c
c########################################################

      character*4 atm1,atm2,atm3
      character*4 res1,res2,res3
      integer i,itmpl
      
      findMolTemplateUreyBradley = 0
      if (NMolTemplateUreyBradleys(itmpl).le.0) return
      
      call sortMolTemplateUreyBradley(atm1,res1,atm2,res2,atm3,res3)
      
      do i=1,NMolTemplateUreyBradleys(itmpl)
        if (atm1.eq.MolTemplateUByName(itmpl,i,1) .and. res1.eq.MolTemplateUByRes(itmpl,i,1)) then
          if (atm2.eq.MolTemplateUByName(itmpl,i,2) .and. res2.eq.MolTemplateUByRes(itmpl,i,2)) then
            if (atm3.eq.MolTemplateUByName(itmpl,i,3) .and. res3.eq.MolTemplateUByRes(itmpl,i,3)) then
              findMolTemplateUreyBradley = i
              return      
            endif
          endif
        endif      
      enddo
      return
      end
