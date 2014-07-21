      double precision function deg2rad(deg)
c################################################
c
c     This is a very simple function that converts
c     deg degrees to radians...
c
c     (C) Christopher Woods 22 January 2002
c
c###############################################

      double precision deg,twopiO360
      parameter(twopiO360=0.017453292519943)
      
      deg2rad = deg * twopiO360
      
      return 
      
      end
