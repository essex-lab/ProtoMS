
c   This file holds the variables associated with pressure (NPT) and volume
      double precision Volume,maxVolChange
      double precision pressure,intPressure
      double precision AvgVolume,AvgVolume2
      integer NVolTry,NVolAcp
      logical nptsim
      
      common /VOLUMEVARS/ Volume,pressure,intPressure,AvgVolume,maxVolChange,AvgVolume2,
     .                    NVolTry,NVolAcp,nptsim

      save /VOLUMEVARS/
      
