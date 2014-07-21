      subroutine printaveragegbenergy
      include 'dimensions.inc'
      include 'constants.inc'
      include 'printer.inc'
      include 'gb.inc'
c###########################################
c     Print running average of the Generalised
c     born energy
c     Julien Michel November 2004
c###########################################
      double precision stdDev
      double precision std,stdf,stdb

      std = stdDev(AvgGBEnergy2,AvgGBEnergy)
      stdf = stdDev(AvgGBEnergyF2,AvgGBEnergyF)
      stdb = stdDev(AvgGBEnergyB2,AvgGBEnergyB)
c      write(printstring,*) 'Average Generalised Born Energy'
c      call prettyLine(RESULTS,printstring)
      write(printstring,45) AvgGBEnergy,std,AvgGBEnergyF,stdf,AvgGBEnergyB,stdb
 45   format('GB Energy    ',f20.10,' ( ',f10.3,' ) | ',f20.10,' ( ',f10.3,' ) ',f20.10,' ( ',f10.3,' )')
      call printLine(RESULTS,printstring)          

      return
      end
