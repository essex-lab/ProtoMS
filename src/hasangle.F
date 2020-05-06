      logical function hasAngle(iprot,zat1,zat2,zat3)
      include 'dimensions.inc'
      include 'proteins.inc'
c###################################################
c
c     This function returns whether iprot contains
c      an angle between zatoms zat1-zat2-zat3
c
c###################################################

      integer i,iprot,zat1,zat2,zat3,ires,nangs,iang

      hasAngle = .false.
      
      ires = Zid2ResID(iprot,zat1,1)
      nangs = NResAngles(iprot,ires)
      if (nangs.le.0) return
      
      hasAngle = .true.
      do i=1,nangs
        iang = ResidueAngle(iprot,ires,i)
        if (Angle(iprot,iang,2).eq.zat2) then
          if (Angle(iprot,iang,1).eq.zat1 .and. Angle(iprot,iang,3).eq.zat3) return
          if (Angle(iprot,iang,1).eq.zat3 .and. Angle(iprot,iang,3).eq.zat1) return
        endif
      enddo

      hasAngle = .false.
      return
      end
