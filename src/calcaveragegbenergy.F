      subroutine calcaveragegbenergy(bigratio,smallratio)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'gb.inc'
c##############################################
c     Running averages of Generalised Born energy
c
c     Julien Michel November 2004
c##############################################
      double precision bigratio,smallratio

      AvgGBEnergy = (AvgGBEnergy*bigratio) + (OldGBEnergy*smallratio)
      AvgGBEnergyF = (AvgGBEnergyF*bigratio) + (OldGBEnergyF*smallratio)
      AvgGBEnergyB = (AvgGBEnergyB*bigratio) + (OldGBEnergyB*smallratio)
      
      AvgGBEnergy2 = (AvgGBEnergy2*bigratio) + (smallratio*(OldGBEnergy**2))
      AvgGBEnergyF2 = (AvgGBEnergyF2*bigratio) + (smallratio*(OldGBEnergyF**2))
      AvgGBEnergyB2 = (AvgGBEnergyB2*bigratio) + (smallratio*(OldGBEnergyB**2))

      return
      end
