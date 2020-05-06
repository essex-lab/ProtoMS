      double precision function ranrange(strtval,endval)
      implicit none
c#####################################################
c
c     This function generates a random dble between
c     strtval and endval
c
c     (C) Christopher Woods   21 August 2002
c
c#####################################################
            
      double precision temp1,temp2
      double precision strtval,endval
      double precision grnd
      
     
      
      if (strtval.gt.endval) then
        temp1 = endval
        temp2 = strtval
      else if (strtval.eq.endval) then
        temp1 = strtval
        temp2 = endval
        ranrange = strtval
        return
      else
        temp1 = strtval
        temp2 = endval
      endif
      
c     use the Mersenne Twister excellent linear random number generator
      ranrange = ( (temp2-temp1) * grnd()) + temp1
      
      return
      end
