       double precision function stdDev(meansq,mean)
       implicit none
c#############################################################
c
c      This calculates the standard deviation, from 
c      the mean of the squares (meansq) and the mean
c
c      (C) Christopher Woods    16 August 2002
c
c##############################################################
       
       double precision meansq,mean

       stdDev = meansq - mean**2
       stdDev = dsqrt(dabs(stdDev))
       
       return
       end
