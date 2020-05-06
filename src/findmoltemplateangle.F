      integer function findMolTemplateAngle(atm1,res1,atm2,res2,atm3,res3,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c########################################################
c
c     this function is the angle equivelent of findMolTemplateBond
c
c     (C) Christopher Woods  18 February 2002
c
c########################################################

      character*4 atm1,atm2,atm3,res1,res2,res3
      integer i,itmpl
      
      findMolTemplateAngle = 0
      if (NMolTemplateAngles(itmpl).le.0) return
      
      call sortMolTemplateAngle(atm1,res1,atm2,res2,atm3,res3)
      
      do i=1,NMolTemplateAngles(itmpl)
        if (atm1.eq.MolTemplateAngName(itmpl,i,1).and.res1.eq.MolTemplateAngRes(itmpl,i,1)) then
          if (atm2.eq.MolTemplateAngName(itmpl,i,2).and.res2.eq.MolTemplateAngRes(itmpl,i,2)) then
            if (atm3.eq.MolTemplateAngName(itmpl,i,3).and.res3.eq.MolTemplateAngRes(itmpl,i,3)) then
              findMolTemplateAngle = i
              return
            endif
          endif
        endif
      enddo
      return
      end
