
      subroutine calcAverageExtraEnergies(bigratio,smallratio)
      include 'dimensions.inc'
      include 'flags.inc'
      include 'constants.inc'
      include 'extraenergies.inc'
      include 'printer.inc'
c#########################################################
c
c     This is the master routine to calculate the average
c     of all of the extra energies. Add your averaging
c     code to this routine...
c
c     (C) Christopher Woods, November 2004
c
c#########################################################

      double precision bigratio,smallratio

c     calculate the average total extra energies      
c      write(printstring,*) 'Avg Extra Nrg was ',avgExtraEnergy, 'bigra' ,bigratio
c     .  ,' smallratio ',smallratio,' old XTRA ',oldExtraEnergy
c      call printLine(DEBUG,printstring)
      avgExtraEnergy = (avgExtraEnergy*bigratio) + (oldExtraEnergy*smallratio)
c      write(printstring,*) 'Avg Extra Nrg is now ',avgExtraEnergy
c      call printLine(DEBUG,printstring)
      avgExtraEnergy2 = (avgExtraEnergy2*bigratio) + (oldExtraEnergy**2 * smallratio)
      avgExtraEnergyB = (avgExtraEnergyB*bigratio) + (oldExtraEnergyB*smallratio)
      avgExtraEnergyB2 = (avgExtraEnergyB2*bigratio) + (oldExtraEnergyB**2 * smallratio)
      avgExtraEnergyF = (avgExtraEnergyF*bigratio) + (oldExtraEnergyF*smallratio)
      avgExtraEnergyF2 = (avgExtraEnergyF2*bigratio) + (oldExtraEnergyF**2 * smallratio)

c     calculate the average solvent cap energies            
      if (SolventCap) then
c        write(printstring,*) 'Avg Svn Cap was ',avgSvnCapEnergy,' Big ',bigratio,
c     .    ' smallratio ',smallratio,' old Cap ',oldSvnCapEnergy
c        call printLine(DEBUG,printstring)
        avgSvnCapEnergy = (avgSvnCapEnergy*bigratio) + (oldSvnCapEnergy * smallratio)
c        write(printstring,*) 'Avg Svn Cap is now ',avgSvnCapEnergy
c        call printLine(DEBUG,printstring)
        avgSvnCapEnergy2 = (avgSvnCapEnergy2*bigratio) + (oldSvnCapEnergy**2 * smallratio)
      endif

c     calculate the average harmonic restraint energies
      if (NHarmonicRestraints.gt.0) then
        avgHarmEnergy = (avgHarmEnergy*bigratio) + (oldHarmEnergy * smallratio)
        avgHarmEnergy2 = (avgHarmEnergy2*bigratio) + (oldHarmEnergy**2 * smallratio)
        avgHarmEnergyB = (avgHarmEnergyB*bigratio) + (oldHarmEnergyB * smallratio)
        avgHarmEnergyB2 = (avgHarmEnergyB2*bigratio) + (oldHarmEnergyB**2 * smallratio)
        avgHarmEnergyF = (avgHarmEnergyF*bigratio) + (oldHarmEnergyF * smallratio)
        avgHarmEnergyF2 = (avgHarmEnergyF2*bigratio) + (oldHarmEnergyF**2 * smallratio)
      endif

      if (NBndHarmonicRestraints.gt.0) then
        avgBndHarmEnergy = (avgBndHarmEnergy*bigratio) + (oldBndHarmEnergy * smallratio)
        avgBndHarmEnergy2 = (avgBndHarmEnergy2*bigratio) + (oldBndHarmEnergy**2 * smallratio)
        avgBndHarmEnergyB = (avgBndHarmEnergyB*bigratio) + (oldBndHarmEnergyB * smallratio)
        avgBndHarmEnergyB2 = (avgBndHarmEnergyB2*bigratio) + (oldBndHarmEnergyB**2 * smallratio)
        avgBndHarmEnergyF = (avgBndHarmEnergyF*bigratio) + (oldBndHarmEnergyF * smallratio)
        avgBndHarmEnergyF2 = (avgBndHarmEnergyF2*bigratio) + (oldBndHarmEnergyF**2 * smallratio)
      endif

c     calculate the average flat bottom restraint energies
      if (NFlatBottomRestraints.gt.0) then
c        write(printstring,*) 'FBOT Add ',oldFBotEnergy,' to ',avgFBotEnergy
c        call printLine(DEBUG,printstring)
        avgFBotEnergy = (avgFBotEnergy*bigratio) + (oldFBotEnergy * smallratio)
        avgFBotEnergy2 = (avgFBotEnergy2*bigratio) + (oldFBotEnergy**2 * smallratio)
        avgFBotEnergyB = (avgFBotEnergyB*bigratio) + (oldFBotEnergyB * smallratio)
        avgFBotEnergyB2 = (avgFBotEnergyB2*bigratio) + (oldFBotEnergyB**2 * smallratio)
        avgFBotEnergyF = (avgFBotEnergyF*bigratio) + (oldFBotEnergyF * smallratio)
        avgFBotEnergyF2 = (avgFBotEnergyF2*bigratio) + (oldFBotEnergyF**2 * smallratio)
      endif
c     calculate the average Dihedral harmonic restraint energies
      if (NDihedralHarmonicRestraints.gt.0) then
        avgDihHarmEnergy = (avgDihHarmEnergy*bigratio) + (oldDihHarmEnergy * smallratio)
        avgDihHarmEnergy2 = (avgDihHarmEnergy2*bigratio) + (oldDihHarmEnergy**2 * smallratio)
        avgDihHarmEnergyB = (avgDihHarmEnergyB*bigratio) + (oldDihHarmEnergyB * smallratio)
        avgDihHarmEnergyB2 = (avgDihHarmEnergyB2*bigratio) + (oldDihHarmEnergyB**2 * smallratio)
        avgDihHarmEnergyF = (avgDihHarmEnergyF*bigratio) + (oldDihHarmEnergyF * smallratio)
        avgDihHarmEnergyF2 = (avgDihHarmEnergyF2*bigratio) + (oldDihHarmEnergyF**2 * smallratio)
      endif   
c     calculate the average Dihedral Lam Torsion restraint energies
      if (NDihedralLamTorRestraints.gt.0) then
        avgDihLamTorEnergy = (avgDihLamTorEnergy*bigratio) + (oldDihLamTorEnergy * smallratio)
        avgDihLamTorEnergy2 = (avgDihLamTorEnergy2*bigratio) + (oldDihLamTorEnergy**2 * smallratio)
        avgDihLamTorEnergyB = (avgDihLamTorEnergyB*bigratio) + (oldDihLamTorEnergyB * smallratio)
        avgDihLamTorEnergyB2 = (avgDihLamTorEnergyB2*bigratio) + (oldDihLamTorEnergyB**2 * smallratio)
        avgDihLamTorEnergyF = (avgDihLamTorEnergyF*bigratio) + (oldDihLamTorEnergyF * smallratio)
        avgDihLamTorEnergyF2 = (avgDihLamTorEnergyF2*bigratio) + (oldDihLamTorEnergyF**2 * smallratio)
      endif    

c     calculate the average hardwall energies
      if (NHardwallPotential.gt.0) then
        avgHardwallEnergy = (avgHardwallEnergy*bigratio) + (oldHardwallEnergy * smallratio)
        avgHardwallEnergy2 = (avgHardwallEnergy2*bigratio) + (oldHardwallEnergy**2 * smallratio)
        avgHardwallEnergyB = (avgHardwallEnergyB*bigratio) + (oldHardwallEnergyB * smallratio)
        avgHardwallEnergyB2 = (avgHardwallEnergyB2*bigratio) + (oldHardwallEnergyB**2 * smallratio)
        avgHardwallEnergyF = (avgHardwallEnergyF*bigratio) + (oldHardwallEnergyF * smallratio)
        avgHardwallEnergyF2 = (avgHardwallEnergyF2*bigratio) + (oldHardwallEnergyF**2 * smallratio)
      endif

c     JM calculate the averages of the surface area energy
      if (doSA) then
        AvgSAEnergy = (AvgSAEnergy*bigratio) + (OldSAEnergy*smallratio)
        AvgSAEnergyF = (AvgSAEnergyF*bigratio) + (OldSAEnergyF*smallratio)
        AvgSAEnergyB = (AvgSAEnergyB*bigratio) + (OldSAEnergyB*smallratio)
        AvgSAEnergy2 = (AvgSAEnergy2*bigratio) + (smallratio*(OldSAEnergy**2))
        AvgSAEnergyF2 = (AvgSAEnergyF2*bigratio) + (smallratio*(OldSAEnergyF**2))
        AvgSAEnergyB2 = (AvgSAEnergyB2*bigratio) + (smallratio*(OldSAEnergyB**2))
        AvgTotalSurface = (AvgTotalSurface*bigratio) + (OldTotalSurface*smallratio)
        AvgTotalSurfaceF = (AvgTotalSurfaceF*bigratio) + (OldTotalSurfaceF*smallratio)
        AvgTotalSurfaceB = (AvgTotalSurfaceB*bigratio) + (OldTotalSurfaceB*smallratio)
        AvgTotalSurface2 = (AvgTotalSurface2*bigratio) + (smallratio*(OldTotalSurface**2))
        AvgTotalSurfaceF2 = (AvgTotalSurfaceF2*bigratio) + (smallratio*(OldTotalSurfaceF**2))
        AvgTotalSurfaceB2 = (AvgTotalSurfaceB2*bigratio) + (smallratio*(OldTotalSurfaceB**2))        
      endif

      if (doGB) then 
c        write(printstring,*) 'Adding to AvgGBEneryg ',OldGBEnergy
c        call printLine(DEBUG,printstring)
        AvgGBEnergy = (AvgGBEnergy*bigratio) + (OldGBEnergy*smallratio)
        AvgGBEnergyF = (AvgGBEnergyF*bigratio) + (OldGBEnergyF*smallratio)
        AvgGBEnergyB = (AvgGBEnergyB*bigratio) + (OldGBEnergyB*smallratio)
        AvgGBEnergy2 = (AvgGBEnergy2*bigratio) + (smallratio*(OldGBEnergy**2))
        AvgGBEnergyF2 = (AvgGBEnergyF2*bigratio) + (smallratio*(OldGBEnergyF**2))
        AvgGBEnergyB2 = (AvgGBEnergyB2*bigratio) + (smallratio*(OldGBEnergyB**2))
      endif

      return
      end
      
