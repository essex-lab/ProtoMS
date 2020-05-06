      subroutine getPerpendicular(xAxis,yAxis)
      include 'dimensions.inc'
      include 'constants.inc'
c#######################################################
c
c     This routine returns in yAxis a vector that
c     is perpendicular to xAxis (it chooses one
c     of the many infinite numbers of vectors!)
c
c     (C) Christopher Woods   4 March 2002
c
c#######################################################

      integer i,iel,iel2
      double precision xAxis(3),yAxis(3),zAxis(3)
      
c     create a copy of xAxis which has two elements swapped
      iel = 0
      iel2 = 0
      do i=3,1,-1
        if (xAxis(i).ne.ZERO) iel = i
        zAxis(i) = xAxis(i)
      enddo

      if (iel.eq.0) then
        do i=1,3
          xAxis(i) = ZERO
          yAxis(i) = ZERO
        enddo
        xAxis(1) = ONE
        yAxis(2) = ONE
        return
      endif
      
      if (iel.eq.1) then
        iel2 = 3
      else if (iel.eq.2) then
        iel2 = 1
      else if (iel.eq.3) then
        iel2 = 2
      endif
      
      zAxis(iel2) = xAxis(iel)
      zAxis(iel) = -xAxis(iel2)
      
c     now construct the yAxis from the cross product of x and z
      call cross(xAxis,zAxis,yAxis)
            
      return
      end
