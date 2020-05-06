      subroutine printEnergies(stp)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'energies.inc'
c#####################################################
c
c     This routine prints out the current energies
c     to the ENERGY stream
c
c     (C) Christopher Woods, October 2004
c
c#####################################################

      integer stp
      
c     print total energies...
      write(printstring,10) stp,oldTotalEnergy,oldTotalEnergyB,oldTotalEnergyF
10    format(i10," Total energy = ",g20.10," backwards = ",g20.10," forwards = ",g20.10)
      call printLine(ENERGY,printstring)            
      
c     need to somehow filter the output...
            
      return
      end
