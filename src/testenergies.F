      subroutine testEnergies
      include 'dimensions.inc'
      include 'enums.inc'
      include 'energies.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'printer.inc'
c#####################################################
c
c     This routine recalculates the total energy
c     of the system from scratch and checks that
c     it agrees with what is currently stored.
c
c     This routine is used to check for energy leaks
c     and can be safely called at any time as it
c     will not write to any global arrays or values
c
c     (C) Christopher Woods, October 2004
c
c#####################################################

      double precision inrg,inrgf,inrgb,inrgf2,inrgb2
      double precision totalEnergy,extraEnergy
      logical ok

c     need to set all moved to calculate total extra energy
      call setAllMoved
            
      inrg = totalEnergy(DONTCARE,inrgb,inrgf) + extraEnergy(DONTCARE,inrgb2,inrgf2)
      inrgb = inrgb + inrgb2
      inrgf = inrgf + inrgf2

c     can now clear the move list
      call clearMoveList      
            
      ok = .true.
      if (dabs(inrg - OldTotalEnergy).gt.SMALL) then
        write(printstring,10) "total energy",inrg,OldTotalEnergy
        call printLine(WARNING,printstring)
        ok = .false.
      endif
10    format("Error in ",a,". It should be ",g20.8," but it is ",g20.8)
            
      if (differentLambdas) then
        if (dabs(inrgf - OldTotalEnergyF).gt.SMALL) then
          write(printstring,10) "total forwards energy",inrg,OldTotalEnergy
          call printLine(WARNING,printstring)
          ok = .false.
        endif

        if (dabs(inrgb - OldTotalEnergyB).gt.SMALL) then
          write(printstring,10) "total backwards energy",inrg,OldTotalEnergy
          call printLine(WARNING,printstring)
          ok = .false.
        endif
      endif
      
      if (ok) return
      
c------------------------------
c    Now down here we should do a more thorough investigation to
c    discover which energy components are wrong...
c     call checkTotalEnergy
      
      return
      end
