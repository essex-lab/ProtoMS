      subroutine getAxisVectors(xAxis,ivec,yAxis,zAxis)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'constants.inc'
c#############################################################
c
c     Given an xaxis, and another vector (ivec) this returns
c     a y-axis that is perpendicular to xAxis, but in the
c     plane of xAxis/ivex, and a zAxis that is perpendicular
c     to both the xAxis and yAxis
c
c     (C) Christopher Woods     19 Febraury 2002
c
c#############################################################

      double precision xAxis(3),ivec(3),yAxis(3),zAxis(3)
      integer i
      double precision dot,veclength2
      double precision alpha,beta
      
c     first check that the vectors are not colinear!
      if (abs(dot(xAxis,ivec)).gt.0.9) then
c       if the axes are near colinear, then it doesn't really matter
c       where the y axis is, as long as it is perpendicular
        
        call getPerpendicular(xAxis,yAxis)
      else
c       copy ivec into yAxis, then normalise it 
        do i=1,3
          yAxis(i) = ivec(i)
        enddo
        call vecnormalise(yAxis)

c       now find the vector perpendicular to 'xAxis' but in plane with xAxis/ivec 
c       and same sign as ivec
c       Get this as  xAxis.yAxis = 0
c           and    alpha*xAxis + beta*yAxis = ivec
c           rearrange to get
c           alpha = (xAxis1.ivec1 + xAxis2.ivec2 + xAxis3.ivec3)
c                                    / Length(xAxis)^2
c       (though Length(xAxis) = 1.0)

        alpha = (xAxis(1)*yAxis(1)+xAxis(2)*yAxis(2)+xAxis(3)*yAxis(3))      
        beta = (veclength2(yAxis) - (alpha**2))
        beta = sqrt(beta)
        if (beta.eq.ZERO) then
          write(printstring,*) "Problem calculating beta in getAxisVectors() beta = 0.0!"
          call printLine(WARNING,printstring)
          call printLine(WARNING,"This is either the result of a programming bug or numerical instability")
          call closeProgram(-1,"ProtoMS is exiting as it may have detected a numerical instability")
        endif
        do i=1,3
          yAxis(i) = (yAxis(i) - (alpha*xAxis(i))) / beta
        enddo
        call vecnormalise(yAxis)

      endif

c     now just get the cross product of the x and y axes
      call cross(xAxis,yAxis,zAxis)
      
      return
      end
