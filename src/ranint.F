      integer function ranInt(strtval,endval)
c#####################################################
c
c     This function generates a random integer between
c     strtval and endval..
c
c     (C) Christopher Woods    4 February 2002
c
c#####################################################

      double precision grnd
      integer strtval,endval,tempstrt,tempend

      ranInt = strtval
      
      if (strtval.eq.endval) return
      if (strtval.gt.endval) then
        tempstrt = endval
        tempend = strtval
      else
        tempstrt = strtval
        tempend = endval
      endif
      
c     generate the random number remembering that 'int' truncates,
c     and thus pretty much zero probability of Ranu = 1.0
c
c     uses Mersenne Twister excellent linear random number generator (c function)
      ranInt = int(grnd()*dble(1+tempend-tempstrt)) + tempstrt
      if (ranInt.gt.tempend) ranInt = tempend
      
      return
      
      end
