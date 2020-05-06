      subroutine updateSolventEnergies
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'energies.inc'
      include 'solutes.inc'
      include 'proteins.inc'
c######################################################
c
c     This routine updates the total energies of the 
c     solvent, after a successful solvent move
c
c     (C) Christopher Woods  11 March 2002
c
c#######################################################

      integer i

      if (NSolvents.le.0) return

c     update the total energy
      oldTotalEnergy = oldTotalEnergy + newIndSvnC + newIndSvnLJ 
     .                                - oldIndSvnC - oldIndSvnLJ

c     update the forwards and backwards energies
      oldTotalEnergyF = oldTotalEnergyF + newIndSvnCF + newIndSvnLJF
     .                                - oldIndSvnCF - oldIndSvnLJF
     
      oldTotalEnergyB = oldTotalEnergyB + newIndSVNCB + newIndSvnLJB
     .                                - oldIndSvnCB - oldIndSvnLJB     
     
      if (NSolvents.gt.1) then     
c       now update the solvent-solvent bulk energy
        OldSvnSvnC = OldSvnSvnC + newIndSvnSvnC - oldIndSvnSvnC
        OldSvnSvnLJ = OldSvnSvnLJ + newIndSvnSvnLJ - oldIndSvnSvnLJ
        OldSvnSvnCF = OldSvnSvnCF + newIndSvnSvnCF - oldIndSvnSvnCF
        OldSvnSvnLJF = OldSvnSvnLJF + newIndSvnSvnLJF - oldIndSvnSvnLJF
        OldSvnSvnCB = OldSvnSvnCB + newIndSvnSvnCB - oldIndSvnSvnCB
        OldSvnSvnLJB = OldSvnSvnLJB + newIndSvnSvnLJB - oldIndSvnSvnLJB
      endif
      
      if (NSolutes.gt.0) then
c       now update the solute-solvent bulk energy
        do 100 i=1,NSolutes
          OldSolSvnC(i) = OldSolSvnC(i) + newIndSvnSolC(i) - oldIndSvnSolC(i)
          OldSolSvnLJ(i) = OldSolSvnLJ(i) + newIndSvnSolLJ(i) - oldIndSvnSolLJ(i)      
          OldSolSvnCF(i) = OldSolSvnCF(i) + newIndSvnSolCF(i) - oldIndSvnSolCF(i)
          OldSolSvnLJF(i) = OldSolSvnLJF(i) + newIndSvnSolLJF(i) - oldIndSvnSolLJF(i)      
          OldSolSvnCB(i) = OldSolSvnCB(i) + newIndSvnSolCB(i) - oldIndSvnSolCB(i)
          OldSolSvnLJB(i) = OldSolSvnLJB(i) + newIndSvnSolLJB(i) - oldIndSvnSolLJB(i)      
100     continue
      endif
      
      if (NProteins.gt.0) then
c       now update the protein-solvent bulk energy
        do 200 i=1,NProteins
          OldProSvnC(i) = OldProSvnC(i) + newIndSvnProC(i) - oldIndSvnProC(i)
          OldProSvnLJ(i) = OldProSvnLJ(i) + newIndSvnProLJ(i) - oldIndSvnProLJ(i)
          OldProSvnCF(i) = OldProSvnCF(i) + newIndSvnProCF(i) - oldIndSvnProCF(i)
          OldProSvnLJF(i) = OldProSvnLJF(i) + newIndSvnProLJF(i) - oldIndSvnProLJF(i)
          OldProSvnCB(i) = OldProSvnCB(i) + newIndSvnProCB(i) - oldIndSvnProCB(i)
          OldProSvnLJB(i) = OldProSvnLJB(i) + newIndSvnProLJB(i) - oldIndSvnProLJB(i)
200     continue
      endif
     
      if (NGCSolutes.gt.0) then
c       now update the GCsolute-solvent bulk energy
        do 300 i=1,NGCSolutes
          OldGCSSvnC(i) = OldGCSSvnC(i) + newIndSvnGCSC(i) - oldIndSvnGCSC(i)
          OldGCSSvnLJ(i) = OldGCSSvnLJ(i) + newIndSvnGCSLJ(i) - oldIndSvnGCSLJ(i)
          OldGCSSvnCF(i) = OldGCSSvnCF(i) + newIndSvnGCSCF(i) - oldIndSvnGCSCF(i)
          OldGCSSvnLJF(i) = OldGCSSvnLJF(i) + newIndSvnGCSLJF(i) - oldIndSvnGCSLJF(i)
          OldGCSSvnCB(i) = OldGCSSvnCB(i) + newIndSvnGCSCB(i) - oldIndSvnGCSCB(i)
          OldGCSSvnLJB(i) = OldGCSSvnLJB(i) + newIndSvnGCSLJB(i) - oldIndSvnGCSLJB(i)
300     continue
      endif
      
      return
      end
