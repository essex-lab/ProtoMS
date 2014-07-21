      subroutine printaveragesaenergy
      include 'dimensions.inc'
      include 'constants.inc'
      include 'printer.inc'
      include 'surface.inc'
c###########################################
c     Print running average of the surface
c     area energy
c     Julien Michel November 2004
c###########################################
      double precision stdDev
      double precision std,stdf,stdb

      std = stdDev(AvgSAEnergy2,AvgSAEnergy)
      stdf = stdDev(AvgSAEnergyF2,AvgSAEnergyF)
      stdb = stdDev(AvgSAEnergyB2,AvgSAEnergyB)
c      write(printstring,*) 'Average Surface Area Energy'
c      call prettyLine(RESULTS,printstring)
      write(printstring,45) AvgSAEnergy,std,AvgSAEnergyF,stdf,AvgSAEnergyB,stdb
 45   format('SA Energy    ',f20.10,' ( ',f10.3,' ) | ',f20.10,' ( ',f10.3,' ) ',f20.10,' ( ',f10.3,' )')
      call printLine(RESULTS,printstring)          

      return
      end
