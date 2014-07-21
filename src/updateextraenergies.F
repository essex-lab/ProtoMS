
      subroutine updateExtraEnergies
      include 'dimensions.inc'
      include 'energies.inc'
      include 'constants.inc'
      include 'extraenergies.inc'
      include 'flags.inc'
      include 'printer.inc'
      include 'quickpot.inc'
c#############################################################################
c
c     This routine update the stored extra energies after
c     a successful move
c
c     (C) Christopher Woods, November 2004
c
c############################################################################

c     update the total system energies      
      oldTotalEnergy = oldTotalEnergy + newExtraEnergyPart - oldExtraEnergyPart
      oldTotalEnergyB = oldTotalEnergyB + newExtraEnergyPartB - oldExtraEnergyPartB
      oldTotalEnergyF = oldTotalEnergyF + newExtraEnergyPartF - oldExtraEnergyPartF
c     update the total extra energies      
c      write(printstring,*) 'Updating XTRA old ',oldExtraEnergy,
c     .  ' + ',newExtraEnergyPart,' - ',oldExtraEnergyPart
c      call printLine(DEBUG,printstring)
      oldExtraEnergy = oldExtraEnergy + newExtraEnergyPart - oldExtraEnergyPart
      oldExtraEnergyB = oldExtraEnergyB + newExtraEnergyPartB - oldExtraEnergyPartB
      oldExtraEnergyF = oldExtraEnergyF + newExtraEnergyPartF - oldExtraEnergyPartF
c      write(printstring,*) 'XTRA now ',oldExtraEnergy
c      call printLine(DEBUG,printstring)

c
c     SA energy 
c
      if (doSA) then
        oldSAEnergy = NewSAEnergy
        oldSAEnergyF = NewSAEnergyF
        oldSAEnergyB = NewSAEnergyB
        oldTotalSurface = NewTotalSurface
        oldTotalSurfaceF = NewTotalSurfaceF
        oldTotalSurfaceB = NewTotalSurfaceB
        call updateSurfaces
      endif
c
c     GB energy
c
      if (doGB) then 
        OldGBEnergy = NewGBEnergy
        OldGBEnergyF = NewGBEnergyF
        OldGBEnergyB = NewGBEnergyB        
        if (.not. ConstBorn) call updateBornRadii
      endif

c     update the solvent cap energies            
      if (SolventCap) then
c        write(printstring,*) 'Updating cap old ',oldSvnCapEnergy,
c     .    ' + ',newSvnCapEnergyPart,' - ',oldSvnCapEnergyPart
c        call printLine(DEBUG,printstring)
        oldSvnCapEnergy = oldSvnCapEnergy + newSvnCapEnergyPart - oldSvnCapEnergyPart
      endif

c     update the harmonic restraint energies
      if (NHarmonicRestraints.gt.0) then
        oldHarmEnergy = oldHarmEnergy + newHarmEnergyPart - oldHarmEnergyPart
        oldHarmEnergyB = oldHarmEnergyB + newHarmEnergyPartB - oldHarmEnergyPartB
        oldHarmEnergyF = oldHarmEnergyF + newHarmEnergyPartF - oldHarmEnergyPartF
      endif  
c     update the bond harmonic restraint energies
      if (NBndHarmonicRestraints.gt.0) then
        oldBndHarmEnergy = oldBndHarmEnergy + newBndHarmEnergyPart - oldBndHarmEnergyPart
        oldBndHarmEnergyB = oldBndHarmEnergyB + newBndHarmEnergyPartB - oldBndHarmEnergyPartB
        oldBndHarmEnergyF = oldBndHarmEnergyF + newBndHarmEnergyPartF - oldBndHarmEnergyPartF
      endif  

c     update the flat bottom restraint energies
      if (NFlatBottomRestraints.gt.0) then
        oldFBotEnergy = oldFBotEnergy + newFBotEnergyPart - oldFBotEnergyPart
        oldFBotEnergyB = oldFBotEnergyB + newFBotEnergyPartB - oldFBotEnergyPartB
        oldFBotEnergyF = oldFBotEnergyF + newFBotEnergyPartF - oldFBotEnergyPartF
      endif  
c     update the dihedral harmonic restraint energies
      if (NDihedralHarmonicRestraints.gt.0) then
        oldDihHarmEnergy = oldDihHarmEnergy + newDihHarmEnergyPart - oldDihHarmEnergyPart
        oldDihHarmEnergyB = oldDihHarmEnergyB + newDihHarmEnergyPartB - oldDihHarmEnergyPartB
        oldDihHarmEnergyF = oldDihHarmEnergyF + newDihHarmEnergyPartF - oldDihHarmEnergyPartF
      endif 
      if (NDihedralLamTorRestraints.gt.0) then
        oldDihLamTorEnergy = oldDihLamTorEnergy + newDihLamTorEnergyPart - oldDihLamTorEnergyPart
        oldDihLamTorEnergyB = oldDihLamTorEnergyB + newDihLamTorEnergyPartB - oldDihLamTorEnergyPartB
        oldDihLamTorEnergyF = oldDihLamTorEnergyF + newDihLamTorEnergyPartF - oldDihLamTorEnergyPartF
      endif 
c     update hardwalls
      if (NHardwallPotential.gt.0) then
        oldHardwallEnergy = oldHardwallEnergy + newHardwallEnergyPart - oldHardwallEnergyPart
        oldHardwallEnergyB = oldHardwallEnergyB + newHardwallEnergyPartB - oldHardwallEnergyPartB
        oldHardwallEnergyF = oldHardwallEnergyF + newHardwallEnergyPartF - oldHardwallEnergyPartF
      endif


       
      return
      end
