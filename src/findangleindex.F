      integer function findAngleIndex(at1,at2,at3)
      include 'dimensions.inc'
      include 'parameters.inc'
c######################################################
c
c     This function returns the index of the angle
c     between amber atoms at1,at2,at3, or returns
c     0 if the angle does not exist
c
c    (C) Christopher Woods
c
c######################################################

      character*2 at1,at2,at3
      integer i

      findAngleIndex = 0
      
      if (NAngParam.gt.0) then
        call sortAngle(at1,at2,at3)
        do i=1,NAngParam
          if (at1.eq.AngName(i,1)) then
            if (at2.eq.AngName(i,2)) then
              if (at3.eq.AngName(i,3)) then
                findAngleIndex = i
                return
              endif
            endif
          endif
        enddo
      endif

      return
      end
