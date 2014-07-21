      logical function hasSoluteAngle(isol,zat1,zat2,zat3)
      include 'dimensions.inc'
      include 'solutes.inc'
c############################################################
c
c     This function returns whether there is an angle between
c     solute atoms zat1-zat2-zat3 in solute isol
c
c     (C) Christopher Woods   19 February 2002
c
c############################################################
      
      integer isol,zat1,zat2,zat3,nangs,i
      
      hasSoluteAngle = .false.
      nangs = NSolAngles(isol)
      if (nangs.le.0) return
      hasSoluteAngle = .true.
      do i=1,nangs
        if (SolAngle(isol,i,1).eq.zat1 .and. SolAngle(isol,i,2).eq.zat2 
     .               .and. SolAngle(isol,i,3).eq.zat3) return
        if (SolAngle(isol,i,1).eq.zat3 .and. SolAngle(isol,i,2).eq.zat2 
     .               .and. SolAngle(isol,i,3).eq.zat1) return
      enddo
      hasSoluteAngle = .false.
      return
      end
