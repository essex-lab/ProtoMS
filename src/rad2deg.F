      double precision function rad2Deg(rad)
      implicit none
      double precision rad,threesixtyOverTwoPI
      parameter(threesixtyOverTwoPi=57.295779513082)
      
      rad2Deg = rad * threesixtyOverTwoPi
      return
      end
