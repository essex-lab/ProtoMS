      subroutine recalculateGeometries
      include 'dimensions.inc'
      include 'solutes.inc'
c################################################################
c
c     This subroutine recalculates the geometries of 
c     solutes that vary with lambda, thus allowing lambda
c     to change throughout the simulation
c
c     (C) Christopher Woods  5 February 2003
c
c################################################################

      integer i
      
      do i=1,NSolutes
        if (.not.IgnoreSolute(i)) then
          if (NSolVariable(i).gt.0) call synchFEPSolutes(i)
        endif
      enddo

      end
