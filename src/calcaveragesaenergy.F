      subroutine calcaveragesaenergy(bigratio,smallratio)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'surface.inc'
c##############################################
c     Running averages of surface area energy
c
c     Julien Michel November 2004
c##############################################
      double precision bigratio,smallratio

      AvgSAEnergy = (AvgSAEnergy*bigratio) + (OldSAEnergy*smallratio)
      AvgSAEnergyF = (AvgSAEnergyF*bigratio) + (OldSAEnergyF*smallratio)
      AvgSAEnergyB = (AvgSAEnergyB*bigratio) + (OldSAEnergyB*smallratio)
      
      AvgSAEnergy2 = (AvgSAEnergy2*bigratio) + (smallratio*(OldSAEnergy**2))
      AvgSAEnergyF2 = (AvgSAEnergyF2*bigratio) + (smallratio*(OldSAEnergyF**2))
      AvgSAEnergyB2 = (AvgSAEnergyB2*bigratio) + (smallratio*(OldSAEnergyB**2))



      return
      end
