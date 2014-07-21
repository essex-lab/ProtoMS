      subroutine updateSoluteEnergies(isol)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'energies.inc'
      include 'solutes.inc'
      include 'proteins.inc'
c############################################################
c
c     This routine updates the total energies of solute
c     isol, and with the rest of the system, after
c     a succesful solute move. The energies are updated
c     as;
c
c     oldSol??? = newSol??? (since they are calculated as 1)
c
c     (C) Christopher Woods   1 March 2002
c
c#############################################################

      integer isol,i
      
      if (isol.le.0 .or. isol.gt.NSolutes) return
      
c     update the total energy
      oldTotalEnergy = oldTotalEnergy + newSoluteEnergy(isol) - oldSoluteEnergy(isol)

c     update the forwards and backwards energies
      oldTotalEnergyF = oldTotalEnergyF + newSoluteEnergyF(isol) - oldSoluteEnergyF(isol)
      oldTotalEnergyB = oldTotalEnergyB + newSoluteEnergyB(isol) - oldSoluteEnergyB(isol)      
                     
c     now copy all of the new energies to old
      oldSoluteEnergy(isol) = newSoluteEnergy(isol)
      oldSoluteEnergyF(isol) = newSoluteEnergyF(isol)
      oldSoluteEnergyB(isol) = newSoluteEnergyB(isol)

c     solute-solvent energy
      if (NSolvents.gt.0) then
        OldSolSvnC(isol) = NewSolSvnC(isol)
        OldSolSvnLJ(isol) = NewSolSvnLJ(isol)
        OldSolSvnCF(isol) = NewSolSvnCF(isol)
        OldSolSvnLJF(isol) = NewSolSvnLJF(isol)
        OldSolSvnCB(isol) = NewSolSvnCB(isol)
        OldSolSvnLJB(isol) = NewSolSvnLJB(isol)
        OldSolSvnDist(isol) = NewSolSvnDist(isol)
      endif

c     solute-GCS energy
      if (NGCSolutes.gt.0) then
        OldSolGCSC(isol) = NewSolGCSC(isol)
        OldSolGCSLJ(isol) = NewSolGCSLJ(isol)
        OldSolGCSCF(isol) = NewSolGCSCF(isol)
        OldSolGCSLJF(isol) = NewSolGCSLJF(isol)
        OldSolGCSCB(isol) = NewSolGCSCB(isol)
        OldSolGCSLJB(isol) = NewSolGCSLJB(isol)
      endif



c     solute protein energies
      if (NProteins.gt.0) then
        do i=1,NProteins
          OldProSolC(i,isol) = NewProSolC(i,isol)
          OldProSolLJ(i,isol) = NewProSolLJ(i,isol)
          OldProSolCF(i,isol) = NewProSolCF(i,isol)
          OldProSolLJF(i,isol) = NewProSolLJF(i,isol)
          OldProSolCB(i,isol) = NewProSolCB(i,isol)
          OldProSolLJB(i,isol) = NewProSolLJB(i,isol)
          OldProSolDist(i,isol) = NewProSolDist(i,isol)
        enddo
      endif

c     solute-solute energies
      if (NSolutes.gt.1) then
        do i=1,NSolutes
          OldSolSolC(isol,i) = NewSolSolC(isol,i)
          OldSolSolLJ(isol,i) = NewSolSolLJ(isol,i)
          OldSolSolC(i,isol) = NewSolSolC(i,isol)
          OldSolSolLJ(i,isol) = NewSolSolLJ(i,isol)
          OldSolSolCF(isol,i) = NewSolSolCF(isol,i)
          OldSolSolLJF(isol,i) = NewSolSolLJF(isol,i)
          OldSolSolCF(i,isol) = NewSolSolCF(i,isol)
          OldSolSolLJF(i,isol) = NewSolSolLJF(i,isol)
          OldSolSolCB(isol,i) = NewSolSolCB(isol,i)
          OldSolSolLJB(isol,i) = NewSolSolLJB(isol,i)
          OldSolSolCB(i,isol) = NewSolSolCB(i,isol)
          OldSolSolLJB(i,isol) = NewSolSolLJB(i,isol)
          OldSolSolDist(isol,i) = NewSolSolDist(isol,i)
          OldSolSolDist(i,isol) = NewSolSolDist(i,isol)
        enddo
      endif

c     now the solute intramolecular interactions
      OldSolBnd(isol) = NewSolBnd(isol)
      OldSolAng(isol) = NewSolAng(isol)
      OldSolDih(isol) = NewSolDih(isol)
      OldSolUBy(isol) = NewSolUBy(isol)
      OldSolNBC(isol) = NewSolNBC(isol)
      OldSolNBLJ(isol) = NewSolNBLJ(isol)
      OldSolBndF(isol) = NewSolBndF(isol)
      OldSolAngF(isol) = NewSolAngF(isol)
      OldSolDihF(isol) = NewSolDihF(isol)
      OldSolUByF(isol) = NewSolUByF(isol)
      OldSolNBCF(isol) = NewSolNBCF(isol)
      OldSolNBLJF(isol) = NewSolNBLJF(isol)
      OldSolBndB(isol) = NewSolBndB(isol)
      OldSolAngB(isol) = NewSolAngB(isol)
      OldSolDihB(isol) = NewSolDihB(isol)
      OldSolUByB(isol) = NewSolUByB(isol)
      OldSolNBCB(isol) = NewSolNBCB(isol)
      OldSolNBLJB(isol) = NewSolNBLJB(isol)
      
      return
      end
