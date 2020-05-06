      integer function findChainUreyBradley(at1,at2,at3,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c#######################################################
c
c      This function returns the index of the Urey-Bradley
c      between atoms at1 and at2 in the chain template
c      itmpl, or 0 if there is no such Urey-Bradley
c
c      (C) Christopher Woods
c
c#######################################################

      integer itmpl,i
      character*4 at1,at2,at3
      
      findChainUreyBradley = 0
      
      if (NChainUreyBradleys(itmpl).le.0) return
      
      call sortTemplateUreyBradley(at1,at2,at3)
      
      do i=1,NChainUreyBradleys(itmpl)
        if (at1.eq.ChainUByName(itmpl,i,1) .and.
     .      at2.eq.ChainUByName(itmpl,i,2) .and.
     .      at3.eq.ChainUByName(itmpl,i,3)) then
          findChainUreyBradley = i
          return
        endif
      enddo
      
      return
      end
