      subroutine updateGCSoluteEnergies(igcs)
      include 'dimensions.inc'
      include 'energies.inc'
      include 'printer.inc'
      include 'flags.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'solvents.inc'
      include 'simparams.inc'
      include 'enums.inc'
c ### delete later below
      include 'parallel.inc'
          
c######################################################
c
c     This routine updates the total energies of the 
c     GC-solute, after a successful GC-solute move
c
c     (C) MSB  17 September 2009
c
c#######################################################

      integer i,igcs

      if (NGCSolutes.le.0) return

c     update the total energy
      oldTotalEnergy = oldTotalEnergy + newIndGCSC + newIndGCSLJ 
     .                                - oldIndGCSC - oldIndGCSLJ

c     update the forwards and backwards energies
      oldTotalEnergyF = oldTotalEnergyF + newIndGCSCF + newIndGCSLJF
     .                                - oldIndGCSCF - oldIndGCSLJF
     
      oldTotalEnergyB = oldTotalEnergyB + newIndGCSCB + newIndGCSLJB
     .                                - oldIndGCSCB - oldIndGCSLJB     
     
      if (NGCSolutes.gt.1) then     
c       now update the GCS-GCS bulk energy
        OldGCSGCSC = OldGCSGCSC + newIndGCSGCSC - OldGCSGCSC
        OldGCSGCSLJ = OldGCSGCSLJ + newIndGCSGCSLJ - OldGCSGCSLJ
        OldGCSGCSCF = OldGCSGCSCF + newIndGCSGCSCF - OldGCSGCSCF
        OldGCSGCSLJF = OldGCSGCSLJF + newIndGCSGCSLJF - OldGCSGCSLJF
        OldGCSGCSCB = OldGCSGCSCB + newIndGCSGCSCB - OldGCSGCSCB
        OldGCSGCSLJB = OldGCSGCSLJB + newIndGCSGCSLJB - OldGCSGCSLJB
      endif
      
      if (NSolvents.gt.0) then
c       now update the GCS-solvent bulk energy
          OldGCSSvnC(igcs) = OldGCSSvnC(igcs) + newIndGCSSvnC(igcs) - oldIndGCSSvnC(igcs)
          OldGCSSvnLJ(igcs) = OldGCSSvnLJ(igcs) + newIndGCSSvnLJ(igcs) - oldIndGCSSvnLJ(igcs)      
          OldGCSSvnCF(igcs) = OldGCSSvnCF(igcs) + newIndGCSSvnCF(igcs) - oldIndGCSSvnCF(igcs)
          OldGCSSvnLJF(igcs) = OldGCSSvnLJF(igcs) + newIndGCSSvnLJF(igcs) - oldIndGCSSvnLJF(igcs)      
          OldGCSSvnCB(igcs) = OldGCSSvnCB(igcs) + newIndGCSSvnCB(igcs) - oldIndGCSSvnCB(igcs)
          OldGCSSvnLJB(igcs) = OldGCSSvnLJB(igcs) + newIndGCSSvnLJB(igcs) - oldIndGCSSvnLJB(igcs)      
100     continue
      endif
      
      if (NProteins.gt.0) then
c       now update the protein-solvent bulk energy
        do 200 i=1,NProteins
          OldProGCSC(i) = OldProGCSC(i) + newIndGCSProC(i) - oldIndGCSProC(i)
          OldProGCSLJ(i) = OldProGCSLJ(i) + newIndGCSProLJ(i) - oldIndGCSProLJ(i)
          OldProGCSCF(i) = OldProGCSCF(i) + newIndGCSProCF(i) - oldIndGCSProCF(i)
          OldProGCSLJF(i) = OldProGCSLJF(i) + newIndGCSProLJF(i) - oldIndGCSProLJF(i)
          OldProGCSCB(i) = OldProGCSCB(i) + newIndGCSProCB(i) - oldIndGCSProCB(i)
          OldProGCSLJB(i) = OldProGCSLJB(i) + newIndGCSProLJB(i) - oldIndGCSProLJB(i)
200     continue
      endif

      if (NSolutes.gt.0) then
c       now update the GCS-solute bulk energy
        if (myid.eq.1) write(6,*) "UPDATE gcse GCSC (before) ",OldSolGCSC(1)
        do 300 i=1,NSolutes
          OldSolGCSC(i) = OldSolGCSC(i) + newIndGCSSolC(i) - oldIndGCSSolC(i)
          OldSolGCSLJ(i) = OldSolGCSLJ(i) + newIndGCSSolLJ(i) - oldIndGCSSolLJ(i)      
          OldSolGCSCF(i) = OldSolGCSCF(i) + newIndGCSSolCF(i) - oldIndGCSSolCF(i)
          OldSolGCSLJF(i) = OldSolGCSLJF(i) + newIndGCSSolLJF(i) - oldIndGCSSolLJF(i)      
          OldSolGCSCB(i) = OldSolGCSCB(i) + newIndGCSSolCB(i) - oldIndGCSSolCB(i)
          OldSolGCSLJB(i) = OldSolGCSLJB(i) + newIndGCSSolLJB(i) - oldIndGCSSolLJB(i)    
        if (myid.eq.1) write(6,*) "UPDATE gcse GCSC (after) ",OldSolGCSC(1)  
300     continue
      endif

      
      return
      end
