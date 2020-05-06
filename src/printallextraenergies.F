
      subroutine printAllExtraEnergies
      include 'dimensions.inc'
      include 'flags.inc'
      include 'constants.inc'
      include 'extraenergies.inc'
      include 'printer.inc'
c##################################################
c
c     Called from printallenergies.F
c
c     Samuel Genheden August 2014
c
c##################################################

      logical printxtra
      
      write(printstring,*) "Average extra energies"
      call prettyLine(RESULTS,printstring)
c     only print out the extra averages if extra energy
c     components are printed out      
      printxtra = .false.
      
      if (SolventCap) then
        printxtra = .true.
c       print out the average solvent cap energies
        write(printstring,20) OldSvnCapEnergy,ZERO
        call printLine(RESULTS,printstring)
20      format("Average solvent cap energy ",g20.10,' ( ',f8.3,' ) ')        
      endif

      if (NHarmonicRestraints.gt.0) then
        printxtra = .true.
c       print out the average harmonic restraint energies
        write(printstring,25) "HARMONIC",OldHarmEnergy,ZERO,
     .                        OldHarmEnergyF,ZERO,
     .                        OldHarmEnergyB,ZERO
        call printLine(RESULTS,printstring)
        write(printstring,26) OldHarmEnergyF-OldHarmEnergy,OldHarmEnergyB-OldHarmEnergy
        call printLine(RESULTS,printstring)
      endif
      
      if (NBndHarmonicRestraints.gt.0) then
        printxtra = .true.
c       print out the average harmonic restraint energies
        write(printstring,25) "BNDHAR",OldBndHarmEnergy,ZERO,
     .                        OldBndHarmEnergyF,ZERO,
     .                        OldBndHarmEnergyB,ZERO
        call printLine(RESULTS,printstring)
        write(printstring,26) OldBndHarmEnergyF-OldBndHarmEnergy,OldBndHarmEnergyB-OldBndHarmEnergy
        call printLine(RESULTS,printstring)
      endif
c############################################
c
c  If you write extra energy terms, then
c  add the 'print' function/statements
c  here!
c      
      if (NFlatBottomRestraints.gt.0) then
        printxtra = .true.
c       print out the average harmonic restraint energies
        write(printstring,25) "FLAT BOT",OldFBotEnergy,ZERO,
     .                        OldFBotEnergyF,ZERO,
     .                        OldFBotEnergyB,ZERO
        call printLine(RESULTS,printstring)
        write(printstring,26) OldFBotEnergyF-OldFBotEnergy,OldFBotEnergyB-OldFBotEnergy
        call printLine(RESULTS,printstring)
      endif

      if (NDihedralHarmonicRestraints.gt.0) then
        printxtra = .true.
c       print out the average harmonic restraint energies
        write(printstring,25) "DIH-HARM",oldDihHarmEnergy,ZERO,
     .                        OldDihHarmEnergyF,ZERO,
     .                        OldDihHarmEnergyB,ZERO
        call printLine(RESULTS,printstring)
        write(printstring,26) OldDihHarmEnergyF-OldDihHarmEnergy,OldDihHarmEnergyB-OldDihHarmEnergy
        call printLine(RESULTS,printstring)
      endif     
      if (NDihedralLamTorRestraints.gt.0) then
        printxtra = .true.
c       print out the average harmonic restraint energies
        write(printstring,25) "DIH-LTOR",OldDihLamTorEnergy,ZERO,
     .                        OldDihLamTorEnergyF,ZERO,
     .                        OldDihLamTorEnergyB,ZERO
        call printLine(RESULTS,printstring)
        write(printstring,26) OldDihLamTorEnergyF-OldDihLamTorEnergy,OldDihLamTorEnergyB-OldDihLamTorEnergy
        call printLine(RESULTS,printstring)
      endif       

      if (NHardwallPotential.gt.0) then
        printxtra = .true.
c       print out the average hardwall energies
        write(printstring,25) "HARDWALL",OldHardwallEnergy,ZERO,
     .                        OldHardwallEnergyF,ZERO,
     .                        OldHardwallEnergyB,ZERO
        call printLine(RESULTS,printstring)
        write(printstring,26) OldHardwallEnergyF-OldHardwallEnergy,OldHardwallEnergyB-OldHardwallEnergy
        call printLine(RESULTS,printstring)
      endif



      if (doSA) then
        printxtra = .true.
c       print out the average surface area
        write(printstring,25) "SURFACE",OldTotalSurface,ZERO,
     .                        OldTotalSurfaceF,ZERO,
     .                        OldTotalSurfaceB,ZERO
        call printLine(RESULTS,printstring)
        write(printstring,26) OldTotalSurfaceF-OldTotalSurface,OldTotalSurfaceB-OldTotalSurface
        call printLine(RESULTS,printstring)
c       print out the average surface area energies
        write(printstring,25) "SA",OldSAEnergy,ZERO,
     .                        OldSAEnergyF,ZERO,
     .                        OldSAEnergyB,ZERO
        call printLine(RESULTS,printstring)
        write(printstring,26) OldSAEnergyF-OldSAEnergy,OldSAEnergyB-OldSAEnergy
        call printLine(RESULTS,printstring)
      endif

      if (doGB) then
        printxtra = .true.
c       print out the average Generalised Born energies
        write(printstring,25) "GB",OldGBEnergy,ZERO,
     .                        OldGBEnergyF,ZERO,
     .                        OldGBEnergyB,ZERO
        call printLine(RESULTS,printstring)
        write(printstring,26) OldGBEnergyF-OldGBEnergy,OldGBEnergyB-OldGBEnergy
        call printLine(RESULTS,printstring)        
      endif
c
c
c############################################     
25      format(a8,1x,g20.10,' ( ',f6.3,' ) | ',g20.10,' ( ',f6.3,' ) ',g20.10,' ( ',f6.3,' )')            
26      format("DIFF     ",34x,g20.10,12x,g20.10)
c     if we have printed something above, then print out the total
c     extra energy
      if (printxtra) then
        call prettyLine(RESULTS,"Average total extra energy")
        write(printstring,10) OldExtraEnergy,ZERO,
     .                        OldExtraEnergyF,ZERO,
     .                        OldExtraEnergyB,ZERO
        call printLine(RESULTS,printstring)
10      format("EXTRA ",g20.10,' ( ',f8.3,' ) | ',g20.10,' ( ',f8.3,' ) ',g20.10,' ( ',f8.3,' )')

        write(printstring,11) OldExtraEnergyF-OldExtraEnergy,OldExtraEnergyB-OldExtraEnergy
        call printLine(RESULTS,printstring)
11      format("DIFF  ",36x,g20.10,14x,g20.10)
      endif
      
      return
      end
      
