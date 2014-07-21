      integer function findMolTemplateDihedral(atm1,res1,atm2,res2,atm3,res3,atm4,res4,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c#########################################################
c
c     This function is the dihedral equilivent of findMolTemplateBond
c
c     (C) Christopher Woods   18 February 2002
c
c#########################################################

      character*4 atm1,atm2,atm3,atm4
      character*4 res1,res2,res3,res4
      integer itmpl,i
      
      findMolTemplateDihedral = 0
      if (NMolTemplateDihedrals(itmpl).le.0) return
      
      call sortMolTemplateDihedral(atm1,res1,atm2,res2,atm3,res3,atm4,res4)
      
      do i=1,NMolTemplateDihedrals(itmpl)
        if (atm1.eq.MolTemplateDihName(itmpl,i,1).and.res1.eq.MolTemplateDihRes(itmpl,i,1)) then
          if (atm2.eq.MolTemplateDihName(itmpl,i,2).and.res2.eq.MolTemplateDihRes(itmpl,i,2)) then
            if (atm3.eq.MolTemplateDihName(itmpl,i,3).and.res3.eq.MolTemplateDihRes(itmpl,i,3)) then
              if (atm4.eq.MolTemplateDihName(itmpl,i,4).and.res4.eq.MolTemplateDihRes(itmpl,i,4)) then
                 findMolTemplateDihedral = i
                 return
              endif
            endif
          endif
        endif      
      enddo
      return
      end
