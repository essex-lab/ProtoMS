      subroutine updateResidueEnergies(iprot)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'energies.inc'
      include 'solutes.inc'
      include 'proteins.inc'
c######################################################
c
c     This routine updates the total energies of
c     protein iprot after a successful residue move.
c
c     The residue move *must* have put the pre-move (old)
c     energies in oldRes???, and the new versions in
c     newRes???. This routine will then calculate the
c     new protein energyies as
c
c     oldPro??? = oldPro??? + newRes??? - oldRes???
c
c     (C) Christopher Woods   7 February 2002
c
c#####################################################

      integer iprot,i
      

      if (iprot.le.0 .or. iprot.gt.NProteins) return
      
c     update the total energy
      oldTotalEnergy = oldTotalEnergy + newResNrg - oldResNrg
c     update the forwards and backwards energies
      oldTotalEnergyF = oldTotalEnergyF + newResNrgF - oldResNrgF
      oldTotalEnergyB = oldTotalEnergyB + newResNrgB - oldResNrgB                   
     
c     update the protein-protein energies
      if (NProteins.gt.1) then
        do i=1,NProteins
          if (i.ne.iprot) then
            OldProProC(iprot,i) = OldProProC(iprot,i) + newResProC(i) - oldResProC(i)
            OldProProC(i,iprot) = OldProProC(iprot,i)
            OldProProCF(iprot,i) = OldProProCF(iprot,i) + newResProCF(i) - oldResProCF(i)
            OldProProCF(i,iprot) = OldProProCF(iprot,i)
            OldProProCB(iprot,i) = OldProProCB(iprot,i) + newResProCB(i) - oldResProCB(i)
            OldProProCB(i,iprot) = OldProProCB(iprot,i)
          
            OldProProLJ(iprot,i) = OldProProLJ(iprot,i) + newResProL(i) - oldResProL(i)
            OldProProLJ(i,iprot) = OldProProLJ(iprot,i)
            OldProProLJF(iprot,i) = OldProProLJF(iprot,i) + newResProLF(i) - oldResProLF(i)
            OldProProLJF(i,iprot) = OldProProLJF(iprot,i)
            OldProProLJB(iprot,i) = OldProProLJB(iprot,i) + newResProLB(i) - oldResProLB(i)
            OldProProLJB(i,iprot) = OldProProLJB(iprot,i)

            if (newResProDist(i).lt.OldProProDist(iprot,i)) then
              OldProProDist(iprot,i) = newResProDist(i)
              OldProProDist(i,iprot) = newResProDist(i)
            endif
          endif
        enddo
      endif
      
c     update the protein-solute energies
      if (NSolutes.gt.0) then
        do i=1,NSolutes
          OldProSolC(iprot,i) = OldProSolC(iprot,i) + newResSolC(i) - oldResSolC(i)
          OldProSolLJ(iprot,i) = OldProSolLJ(iprot,i) + newResSolL(i) - oldResSolL(i)
          OldProSolCF(iprot,i) = OldProSolCF(iprot,i) + newResSolCF(i) - oldResSolCF(i)
          OldProSolLJF(iprot,i) = OldProSolLJF(iprot,i) + newResSolLF(i) - oldResSolLF(i)
          OldProSolCB(iprot,i) = OldProSolCB(iprot,i) + newResSolCB(i) - oldResSolCB(i)
          OldProSolLJB(iprot,i) = OldProSolLJB(iprot,i) + newResSolLB(i) - oldResSolLB(i)
          if (newResSolDist(i).lt.OldProSolDist(iprot,i))
     .                           OldProSolDist(iprot,i) = newResSolDist(i)
        enddo
      endif
      
c     update the protein-solvent energies
      if (NSolvents.gt.0) then
        OldProSvnC(iprot) = OldProSvnC(iprot) + newResSvnC - oldResSvnC
        OldProSvnLJ(iprot) = OldProSvnLJ(iprot) + newResSvnL - oldResSvnL
        OldProSvnCF(iprot) = OldProSvnCF(iprot) + newResSvnCF - oldResSvnCF
        OldProSvnLJF(iprot) = OldProSvnLJF(iprot) + newResSvnLF - oldResSvnLF
        OldProSvnCB(iprot) = OldProSvnCB(iprot) + newResSvnCB - oldResSvnCB
        OldProSvnLJB(iprot) = OldProSvnLJB(iprot) + newResSvnLB - oldResSvnLB
        if (newResSvnDist.lt.OldProSvnDist(iprot)) OldProSvnDist(iprot) = newResSvnDist
      endif

c     update the protein-GCsolute energies
      if (NGCSolutes.gt.0) then
        OldProGCSC(iprot) = OldProGCSC(iprot) + newResGCSC - oldResGCSC
        OldProGCSLJ(iprot) = OldProGCSLJ(iprot) + newResGCSL - oldResGCSL
        OldProGCSCF(iprot) = OldProGCSCF(iprot) + newResGCSCF - oldResGCSCF
        OldProGCSLJF(iprot) = OldProGCSLJF(iprot) + newResGCSLF - oldResGCSLF
        OldProGCSCB(iprot) = OldProGCSCB(iprot) + newResGCSCB - oldResGCSCB
        OldProGCSLJB(iprot) = OldProGCSLJB(iprot) + newResGCSLB - oldResGCSLB
        if (newResGCSDist.lt.OldProGCSDist(iprot)) OldProGCSDist(iprot) = newResGCSDist
      endif
      
c     update the protein intramolecular terms
      OldProBnd(iprot) = OldProBnd(iprot) + newResBnd - oldResBnd
      OldProAng(iprot) = OldProAng(iprot) + newResAng - oldResAng
      OldProDih(iprot) = OldProDih(iprot) + newResDih - oldResDih
      OldProUBy(iprot) = OldProUBy(iprot) + newResUBy - oldResUBy
      OldProNBC(iprot) = OldProNBC(iprot) + newResNBC - oldResNBC
      OldProNBLJ(iprot) = OldProNBLJ(iprot) + newResNBLJ - oldResNBLJ
      
      OldProBndF(iprot) = OldProBndF(iprot) + newResBndF - oldResBndF
      OldProBndB(iprot) = OldProBndB(iprot) + newResBndB - oldResBndB
      OldProAngF(iprot) = OldProAngF(iprot) + newResAngF - oldResAngF
      OldProAngB(iprot) = OldProAngB(iprot) + newResAngB - oldResAngB
      OldProDihF(iprot) = OldProDihF(iprot) + newResDihF - oldResDihF
      OldProDihB(iprot) = OldProDihB(iprot) + newResDihB - oldResDihB
      OldProUByF(iprot) = OldProUByF(iprot) + newResUByF - oldResUByF
      OldProUByB(iprot) = OldProUByB(iprot) + newResUByB - oldResUByB
      OldProNBCF(iprot) = OldProNBCF(iprot) + newResNBCF - oldResNBCF
      OldProNBCB(iprot) = OldProNBCB(iprot) + newResNBCB - oldResNBCB
      OldProNBLJF(iprot) = OldProNBLJF(iprot) + newResNBLJF - oldResNBLJF
      OldProNBLJB(iprot) = OldProNBLJB(iprot) + newResNBLJB - oldResNBLJB
      
      return
      end
