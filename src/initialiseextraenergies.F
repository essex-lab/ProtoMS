
      subroutine initialiseExtraEnergies(energiesonly)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'extraenergies.inc'

c#########################################################
c
c     This routine is used to initialise all of the 
c     extra energy variables. If 'energiesonly' is 
c     true then you should only initialise the 
c     old/new energy terms, not the avg terms or
c     any other extra energy variables.
c
c     (C) Christopher Woods, November 2004
c
c#########################################################

      logical energiesonly

c     initialise the total extra energies
      oldExtraEnergy = ZERO
      oldExtraEnergyB = ZERO
      oldExtraEnergyF = ZERO
      oldExtraEnergyPart = ZERO
      newExtraEnergyPart = ZERO
      oldExtraEnergyPartB = ZERO
      newExtraEnergyPartB = ZERO
      oldExtraEnergyPartF = ZERO
      newExtraEnergyPartF = ZERO

      if (.not.energiesonly) then
        avgExtraEnergy = ZERO
        avgExtraEnergy2 = ZERO
        avgExtraEnergyB = ZERO
        avgExtraEnergyB2 = ZERO
        avgExtraEnergyF = ZERO
        avgExtraEnergyF2 = ZERO
      endif 
      
c     initialise solvent cap energies
      oldSvnCapEnergy = ZERO
      oldSvnCapEnergyPart = ZERO
      newSvnCapEnergyPart = ZERO

      if (.not.energiesonly) then
        avgSvnCapEnergy = ZERO
        avgSvnCapEnergy2 = ZERO
      endif   
      
c     initialise the harmonic restraint energies
      oldHarmEnergy = ZERO
      oldHarmEnergyB = ZERO
      oldHarmEnergyF = ZERO
      oldHarmEnergyPart = ZERO
      newHarmEnergyPart = ZERO
      oldHarmEnergyPartB = ZERO
      newHarmEnergyPartB = ZERO
      oldHarmEnergyPartF = ZERO
      newHarmEnergyPartF = ZERO
      
      if (.not.energiesonly) then
        NHarmonicRestraints = 0
        avgHarmEnergy = ZERO
        avgHarmEnergy2 = ZERO
        avgHarmEnergyB = ZERO
        avgHarmEnergyB2 = ZERO
        avgHarmEnergyF = ZERO
        avgHarmEnergyF2 = ZERO
      endif  
c     initialise the flat bottom restraint energies
      oldFBotEnergy = ZERO
      oldFBotEnergyB = ZERO
      oldFBotEnergyF = ZERO
      oldFBotEnergyPart = ZERO
      newFBotEnergyPart = ZERO
      oldFBotEnergyPartB = ZERO
      newFBotEnergyPartB = ZERO
      oldFBotEnergyPartF = ZERO
      newFBotEnergyPartF = ZERO
      
      if (.not.energiesonly) then
        NFlatBottomRestraints = 0
        avgFBotEnergy = ZERO
        avgFBotEnergy2 = ZERO
        avgFBotEnergyB = ZERO
        avgFBotEnergyB2 = ZERO
        avgFBotEnergyF = ZERO
        avgFBotEnergyF2 = ZERO
      endif 

c     initialise dihedral harmonic restraints
      oldDihHarmEnergy = ZERO
      oldDihHarmEnergyB = ZERO
      oldDihHarmEnergyF = ZERO
      oldDihHarmEnergyPart = ZERO
      newDihHarmEnergyPart = ZERO
      oldDihHarmEnergyPartB = ZERO
      newDihHarmEnergyPartB = ZERO
      oldDihHarmEnergyPartF = ZERO
      newDihHarmEnergyPartF = ZERO
      
      if (.not.energiesonly) then
        NDihedralHarmonicRestraints = 0
        avgDihHarmEnergy = ZERO
        avgDihHarmEnergy2 = ZERO
        avgDihHarmEnergyB = ZERO
        avgDihHarmEnergyB2 = ZERO
        avgDihHarmEnergyF = ZERO
        avgDihHarmEnergyF2 = ZERO
      endif  
c     initialise dihedral lamtorsion restraints
      oldDihLamTorEnergy = ZERO
      oldDihLamTorEnergyB = ZERO
      oldDihLamTorEnergyF = ZERO
      oldDihLamTorEnergyPart = ZERO
      newDihLamTorEnergyPart = ZERO
      oldDihLamTorEnergyPartB = ZERO
      newDihLamTorEnergyPartB = ZERO
      oldDihLamTorEnergyPartF = ZERO
      newDihLamTorEnergyPartF = ZERO
      
      if (.not.energiesonly) then
        nDihedralLamTorRestraints = 0
        avgDihLamTorEnergy = ZERO
        avgDihLamTorEnergy2 = ZERO
        avgDihLamTorEnergyB = ZERO
        avgDihLamTorEnergyB2 = ZERO
        avgDihLamTorEnergyF = ZERO
        avgDihLamTorEnergyF2 = ZERO
      endif 

c     initialise the hardwall energies
      oldHardwallEnergy = ZERO
      oldHardwallEnergyB = ZERO
      oldHardwallEnergyF = ZERO
      oldHardwallEnergyPart = ZERO
      newHardwallEnergyPart = ZERO
      oldHardwallEnergyPartB = ZERO
      newHardwallEnergyPartB = ZERO
      oldHardwallEnergyPartF = ZERO
      newHardwallEnergyPartF = ZERO

      if (.not.energiesonly) then
        NHardwallPotential = 0
        avgHardwallEnergy = ZERO
        avgHardwallEnergy2 = ZERO
        avgHardwallEnergyB = ZERO
        avgHardwallEnergyB2 = ZERO
        avgHardwallEnergyF = ZERO
        avgHardwallEnergyF2 = ZERO
      endif


c     JM initialise the SA energy
      oldTotalSurface = ZERO
      oldTotalSurfaceF = ZERO
      oldTotalSurfaceB = ZERO
      oldSAEnergy = ZERO
      oldSAEnergyF = ZERO
      oldSAEnergyB = ZERO
      NewTotalSurface = ZERO
      NewTotalSurfaceF = ZERO
      NewTotalSurfaceB = ZERO
      NewSAEnergy = ZERO
      NewSAEnergyF = ZERO
      NewSAEnergyB = ZERO
      if (doSA) then 
        call clearSurfaces
      endif

      if (.not.energiesonly) then
        AvgSAEnergy = ZERO
        AvgSAEnergyF = ZERO
        AvgSAEnergyB = ZERO
        AvgSAEnergy2 = ZERO
        AvgSAEnergyF2 = ZERO
        AvgSAEnergyB2 = ZERO
      endif

c     JM initialise the GB energy
      OldGBEnergy = ZERO
      OldGBEnergyF = ZERO
      OldGBEnergyB = ZERO
      NewGBEnergy = ZERO
      NewGBEnergyF = ZERO
      NewGBEnergyB = ZERO
      if (doGB) then 
        call ZeroAllPDACont
        call ZeroPDACont
        call gbsavesystem
      endif
      if (.not. energiesonly) then
        AvgGBEnergy = ZERO
        AvgGBEnergyF = ZERO
        AvgGBEnergyB = ZERO
        AvgGBEnergy2 = ZERO
        AvgGBEnergyF2 = ZERO
        AvgGBEnergyB2 = ZERO
      endif
      
      return
      end
      
