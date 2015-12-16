      subroutine checkGotTotal
      include 'dimensions.inc'
      include 'printer.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'gb.inc'
c#########################################################
c
c       This routine checks to see if the initial total
c       energy has been calculated, and if not, then 
c       it calculates it!
c
c       (C) Christopher Woods    9 August 2002
c
c#########################################################
      
      double precision nrg,nrgf,nrgb,inrgb,inrgf
      double precision totalEnergy,extraEnergy

      if (GotInitialTotalEnergy) return
                       
      call printLine(INFO,"Calculating the initial total energy from scratch...")

c     need to set all moved to calculate the total extra energies
      call setAllMoved      

c     need to clear the extra energies as they are being recalculated from scratch
      call initialiseExtraEnergies(.true.)      

c     calculate the total energy - note that the 'NEW' extra energy is calculated
c     as this will be used to update the OLD extra energies                        
      nrg = totalEnergy(OLD,nrgb,nrgf) + extraEnergy(NEW,inrgb,inrgf)
      nrgb = nrgb + inrgb
      nrgf = nrgf + inrgf
      
c     clear the move list
      call clearMoveList
      
c     need to update the extra energies to update 'OLD' from 'NEW'
      call updateExtraEnergies
      
      write(printstring,10) nrg,nrgf,nrgb
10    format('Total Energy = ',f15.4,' (',f15.4,' | ',f15.4,')')
      call printLine(INFO,printstring)

      GotInitialTotalEnergy = .true.
      
      return
      end
