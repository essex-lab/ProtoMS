      integer function findBondIndex(at1,at2)
      include 'dimensions.inc'
      include 'parameters.inc'
c######################################################
c
c     This function finds the index of the bond between
c     at1 and at2 in the bond parameter database or 
c     0 if the bond does not exist
c
c     (C) Christopher Woods
c
c######################################################
      
      character*2 at1,at2
      integer i
      
      findBondIndex = 0

      call sortBond(at1,at2)

      if (NBndParam.gt.0) then      
        do i=1,NBndParam
          if (at1.eq.BndName(i,1)) then
            if (at2.eq.BndName(i,2)) then
              findBondIndex = i
              return
            endif
          endif
        enddo
      endif
            
      return
      end
