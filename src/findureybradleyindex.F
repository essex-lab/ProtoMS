      integer function findUreyBradleyIndex(at1,at2,at3)
      include 'dimensions.inc'
      include 'parameters.inc'
c#########################################################
c
c     Find the index of the ub term between atoms
c     with amber types at1 and at2, or 0 if no
c     such parameter exists
c 
c     (C) Christopher Woods
c
c#########################################################

      character*2 at1,at2,at3
      integer i
      
      findUreyBradleyIndex = 0
      
      call sortUreyBradley(at1,at2,at3)

      if (NUByParam.gt.0) then      
        do i=1,NUByParam
          if (at1.eq.UByName(i,1)) then
            if (at2.eq.UByName(i,2)) then
              if (at3.eq.UByName(i,3)) then
                findUreyBradleyIndex = i
                return
              endif
            endif
          endif
        enddo
      endif

      return
      end
