      integer function findMolTemplateBond(atm1,res1,atm2,res2,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c########################################################
c
c     This finds the bonds between at1(res1)-at2(res2) in
c     molecule template itmpl, and returns the index number,
c     or 0 if it does not exist
c
c     (C) Christopher Woods    18 February 2002
c
c########################################################

      character*4 atm1,res1,atm2,res2
      integer itmpl,i
      
      findMolTemplateBond = 0
      if (NMolTemplateBonds(itmpl).le.0) return
      
      call sortMolTemplateBond(atm1,res1,atm2,res2)
      
      do i=1,NMolTemplateBonds(itmpl)
        if (atm1.eq.MolTemplateBndName(itmpl,i,1) .and. res1.eq.MolTemplateBndRes(itmpl,i,1)) then
          if (atm2.eq.MolTemplateBndName(itmpl,i,2) .and. res2.eq.MolTemplateBndRes(itmpl,i,2)) then
            findMolTemplateBond = i
            return      
          endif
        endif      
      enddo
      return
      end
