      integer function findChainBond(at1,at2,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c#######################################################
c
c      This function returns the index of the bond
c      between atoms at1 and at2 in the chain template
c      itmpl, or 0 if there is no such bond
c
c      (C) Christopher Woods
c
c#######################################################

      integer itmpl,i
      character*4 at1,at2
      
      findChainBond = 0
      
      if (NChainBonds(itmpl).le.0) return
      
      call sortTemplateBond(at1,at2)
      
      do i=1,NChainBonds(itmpl)
        if (at1.eq.ChainBndName(itmpl,i,1) .and.
     .      at2.eq.ChainBndName(itmpl,i,2)) then
          findChainBond = i
          return
        endif
      enddo
      
      return
      end
