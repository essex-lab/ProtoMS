
      subroutine printAverageExtraEnergies
      include 'dimensions.inc'
      include 'flags.inc'
      include 'constants.inc'
      include 'extraenergies.inc'
      include 'printer.inc'
c##################################################
c
c     This routine prints out the average 
c     extra energies
c
c     (C) Christopher Woods, November 2004
c
c##################################################

      double precision stdDev
      logical printxtra
      
      write(printstring,*) "Average extra energies"
      call prettyLine(RESULTS,printstring)
c     only print out the extra averages if extra energy
c     components are printed out      
      printxtra = .false.
      
      if (SolventCap) then
        printxtra = .true.
c       print out the average solvent cap energies
        write(printstring,20) avgSvnCapEnergy,stdDev(avgSvnCapEnergy2,avgSvnCapEnergy)
        call printLine(RESULTS,printstring)
20      format("Average solvent cap energy ",g20.10,' ( ',f8.3,' ) ')        
      endif

      if (NHarmonicRestraints.gt.0) then
        printxtra = .true.
c       print out the average harmonic restraint energies
        write(printstring,25) "HARMONIC",avgHarmEnergy,stdDev(avgHarmEnergy2,avgHarmEnergy),
     .                        avgHarmEnergyF,stdDev(avgHarmEnergyF2,avgHarmEnergyF),
     .                        avgHarmEnergyB,stdDev(avgHarmEnergyB2,avgHarmEnergyB)
        call printLine(RESULTS,printstring)
        write(printstring,26) avgHarmEnergyF-avgHarmEnergy,avgHarmEnergyB-avgHarmEnergy
        call printLine(RESULTS,printstring)
      endif
      
      if (NBndHarmonicRestraints.gt.0) then
        printxtra = .true.
c       print out the average harmonic restraint energies
        write(printstring,25) "BNDHAR",avgBndHarmEnergy,stdDev(avgBndHarmEnergy2,avgBndHarmEnergy),
     .                        avgBndHarmEnergyF,stdDev(avgBndHarmEnergyF2,avgBndHarmEnergyF),
     .                        avgBndHarmEnergyB,stdDev(avgBndHarmEnergyB2,avgBndHarmEnergyB)
        call printLine(RESULTS,printstring)
        write(printstring,26) avgBndHarmEnergyF-avgBndHarmEnergy,avgBndHarmEnergyB-avgBndHarmEnergy
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
        write(printstring,25) "FLAT BOT",avgFBotEnergy,stdDev(avgFBotEnergy2,avgFBotEnergy),
     .                        avgFBotEnergyF,stdDev(avgFBotEnergyF2,avgFBotEnergyF),
     .                        avgFBotEnergyB,stdDev(avgFBotEnergyB2,avgFBotEnergyB)
        call printLine(RESULTS,printstring)
        write(printstring,26) avgFBotEnergyF-avgFBotEnergy,avgFBotEnergyB-avgFBotEnergy
        call printLine(RESULTS,printstring)
      endif

      if (NDihedralHarmonicRestraints.gt.0) then
        printxtra = .true.
c       print out the average harmonic restraint energies
        write(printstring,25) "DIH-HARM",avgDihHarmEnergy,stdDev(avgDihHarmEnergy2,avgDihHarmEnergy),
     .                        avgDihHarmEnergyF,stdDev(avgDihHarmEnergyF2,avgDihHarmEnergyF),
     .                        avgDihHarmEnergyB,stdDev(avgDihHarmEnergyB2,avgDihHarmEnergyB)
        call printLine(RESULTS,printstring)
        write(printstring,26) avgDihHarmEnergyF-avgDihHarmEnergy,avgDihHarmEnergyB-avgDihHarmEnergy
        call printLine(RESULTS,printstring)
      endif     
      if (NDihedralLamTorRestraints.gt.0) then
        printxtra = .true.
c       print out the average harmonic restraint energies
        write(printstring,25) "DIH-LTOR",avgDihLamTorEnergy,stdDev(avgDihLamTorEnergy2,avgDihLamTorEnergy),
     .                        avgDihLamTorEnergyF,stdDev(avgDihLamTorEnergyF2,avgDihLamTorEnergyF),
     .                        avgDihLamTorEnergyB,stdDev(avgDihLamTorEnergyB2,avgDihLamTorEnergyB)
        call printLine(RESULTS,printstring)
        write(printstring,26) avgDihLamTorEnergyF-avgDihLamTorEnergy,avgDihLamTorEnergyB-avgDihLamTorEnergy
        call printLine(RESULTS,printstring)
      endif       

      if (NHardwallPotential.gt.0) then
        printxtra = .true.
c       print out the average hardwall energies
        write(printstring,25) "HARDWALL",avgHardwallEnergy,stdDev(avgHardwallEnergy2,avgHardwallEnergy),
     .                        avgHardwallEnergyF,stdDev(avgHardwallEnergyF2,avgHardwallEnergyF),
     .                        avgHardwallEnergyB,stdDev(avgHardwallEnergyB2,avgHardwallEnergyB)
        call printLine(RESULTS,printstring)
        write(printstring,26) avgHardwallEnergyF-avgHardwallEnergy,avgHardwallEnergyB-avgHardwallEnergy
        call printLine(RESULTS,printstring)
      endif



      if (doSA) then
        printxtra = .true.
c       print out the average surface area
        write(printstring,25) "SURFACE",AvgTotalSurface,stdDev(AvgTotalSurface2,AvgTotalSurface),
     .                        AvgTotalSurfaceF,stdDev(AvgTotalSurfaceF2,AvgTotalSurfaceF),
     .                        AvgTotalSurfaceB,stdDev(AvgTotalSurfaceB2,AvgTotalSurfaceB)
        call printLine(RESULTS,printstring)
        write(printstring,26) AvgTotalSurfaceF-AvgTotalSurface,AvgTotalSurfaceB-AvgTotalSurface
        call printLine(RESULTS,printstring)
c       print out the average surface area energies
        write(printstring,25) "SA",AvgSAEnergy,stdDev(AvgSAEnergy2,AvgSAEnergy),
     .                        AvgSAEnergyF,stdDev(AvgSAEnergyF2,AvgSAEnergyF),
     .                        AvgSAEnergyB,stdDev(AvgSAEnergyB2,AvgSAEnergyB)
        call printLine(RESULTS,printstring)
        write(printstring,26) AvgSAEnergyF-AvgSAEnergy,AvgSAEnergyB-AvgSAEnergy
        call printLine(RESULTS,printstring)
      endif

      if (doGB) then
        printxtra = .true.
c       print out the average Generalised Born energies
        write(printstring,25) "GB",AvgGBEnergy,stdDev(AvgGBEnergy2,AvgGBEnergy),
     .                        AvgGBEnergyF,stdDev(AvgGBEnergyF2,AvgGBEnergyF),
     .                        AvgGBEnergyB,stdDev(AvgGBEnergyB2,AvgGBEnergyB)
        call printLine(RESULTS,printstring)
        write(printstring,26) AvgGBEnergyF-AvgGBEnergy,AvgGBEnergyB-AvgGBEnergy
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
        write(printstring,10) avgExtraEnergy,stdDev(avgExtraEnergy2,avgExtraEnergy),
     .                        avgExtraEnergyF,stdDev(avgExtraEnergyF2,avgExtraEnergyF),
     .                        avgExtraEnergyB,stdDev(avgExtraEnergyB2,avgExtraEnergyB)
        call printLine(RESULTS,printstring)
10      format("EXTRA ",g20.10,' ( ',f8.3,' ) | ',g20.10,' ( ',f8.3,' ) ',g20.10,' ( ',f8.3,' )')

        write(printstring,11) avgExtraEnergyF-avgExtraEnergy,avgExtraEnergyB-avgExtraEnergy
        call printLine(RESULTS,printstring)
11      format("DIFF  ",36x,g20.10,14x,g20.10)
      endif
      
      return
      end
      
