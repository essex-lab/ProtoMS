      integer function findChainAngle(at1,at2,at3,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c#######################################################
c
c      This function returns the index of the angle
c      between atoms at1-at2-at3 in the chain template
c      itmpl, or 0 if there is no such angle
c
c      (C) Christopher Woods
c
c#######################################################

      integer itmpl,i
      character*4 at1,at2,at3
      
      findChainAngle = 0
      
      if (NChainAngles(itmpl).le.0) return
      
      call sortTemplateAngle(at1,at2,at3)
      
      do i=1,NChainAngles(itmpl)
        if (at1.eq.ChainAngName(itmpl,i,1) .and.
     .      at2.eq.ChainAngName(itmpl,i,2) .and.
     .      at3.eq.ChainAngName(itmpl,i,3)) then
          findChainAngle = i
          return
        endif
      enddo
      
      return
      end
