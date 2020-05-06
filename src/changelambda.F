      subroutine changeLambda
      include 'dimensions.inc'
      include 'enums.inc'
      include 'flags.inc'
      include 'simparams.inc'
      include 'energies.inc'
      include 'printer.inc'
c####################################################
c
c     This routine calculates the energy for the
c     current value of lambda, printing the change
c     in energy from the old configuration.
c
c     (C) Christopher Woods, October 2004
c
c####################################################

      double precision ieold,ieoldb,ieoldf
      double precision inrg,inrgb,inrgf,inb,inf
      double precision totalEnergy,extraEnergy
      
c     recalculate the geometries and parameters
c     for this value of lambda
      call recalculateParameters()
      call recalculateGeometries()
      
c     now save the old total energy
      ieold = oldTotalEnergy
      ieoldf = oldTotalEnergyF
      ieoldb = oldTotalEnergyB
      
c     see if the lambdas are different
      call scanForVariableGroups()

c     set everything as having moved
      call setAllMoved()
c     clear the extra energies
      call initialiseExtraEnergies(.true.)
      
c     now calculate the total energy 
      inrg = totalEnergy(OLD,inrgb,inrgf) + extraEnergy(NEW,inb,inf)
      inrgb = inrgb + inb
      inrgf = inrgf + inf
      
c     clear the move list and ensure that the extra energies are stored in 'OLD'
      call clearMoveList()
      call updateExtraEnergies()
      
      ieold = inrg - ieold
      
      if (differentLambdas) then
        ieoldb = inrgb - ieoldb
        ieoldf = inrgf - ieoldf
        write(printstring,10) lambda,inrg,ieold,lambdab,ieoldb,lambdaf,ieoldf
10      format("Changed lambda to ",f8.3,". Energy = ",g20.8," Change in energy = ",g20.8,
     .         " (lambdab, ",f8.3," dE = ",g20.8," lambdaf, ",f8.3," dE = ",g20.8,")")
      else
        write(printstring,11) lambda,inrg,ieold
11      format("Changed lambda to ",f8.3,". Energy = ",g20.8," Change in energy = ",g20.8)
      endif
      
      call printLine(INFO,printstring)
      call printLine(ENERGY,printstring)
      
      return
      end
