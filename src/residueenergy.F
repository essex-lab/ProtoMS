      double precision function residueEnergy(savenrg,iprot,ires,nrgb,nrgf)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'simparams.inc'
      include 'enums.inc'
      include 'energies.inc'
c##################################################################
c
c     This function calculates and returns the total energy of 
c     residue ires in protein iprot...
c
c     If savenrg=1 then we are calculating an old energy, if
c     savenrg=2 then this is a new energy. Otherwise, dont save...
c
c     (C) Christopher Woods    4 February 2002
c
c     Modified by Samuel Genheden 2012-2013
c
c##################################################################

      integer ires,iprot,i,savenrg
      double precision bndnrg,angnrg,dihnrg,ubynrg,nbcnrg,nbljnrg
      double precision cnrg,ljnrg,icnrg,iljnrg,mindist
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb,icnrgf,icnrgb,iljnrgf,iljnrgb
      double precision bndnrgf,angnrgf,dihnrgf,ubynrgf,nbcnrgf,nbljnrgf
      double precision bndnrgb,angnrgb,dihnrgb,ubynrgb,nbcnrgb,nbljnrgb
      double precision scl,sclf,sclb,nrg
      double precision getLowScale,getHighScale
      double precision residueIntraEnergy,residueProteinEnergy
      double precision residueSoluteEnergy,residueSolventEnergy,residueGCSoluteEnergy
      double precision nrgf,nrgb

      residueEnergy = ZERO
      cnrgf = ZERO
      ljnrgf = ZERO
      cnrgb = ZERO
      ljnrgb = ZERO
      nrgf = ZERO
      nrgb = ZERO

      if (savenrg.eq.OLD) then
        oldResSvnC = ZERO
        oldResSvnL = ZERO
        oldResSvnCF = ZERO
        oldResSvnLF = ZERO
        oldResSvnCB = ZERO
        oldResSvnLB = ZERO
        oldResGCSC = ZERO
        oldResGCSL = ZERO
        oldResGCSCF = ZERO
        oldResGCSLF = ZERO
        oldResGCSCB = ZERO
        oldResGCSLB = ZERO
        oldResBnd = ZERO
        oldResAng = ZERO
        oldResDih = ZERO
        oldResUBy = ZERO
        oldResNBC = ZERO
        oldResNBLJ = ZERO
        oldResSvnDist = ZERO
        oldResGCSDist = ZERO
        oldResNrg = ZERO 
        oldResNrgF = ZERO
        oldResNrgB = ZERO
        oldResBndF = ZERO
        oldResBndB = ZERO
        oldResAngF = ZERO
        oldResAngB = ZERO
        oldResDihF = ZERO
        oldResDihB = ZERO
        oldResUByF = ZERO
        oldResUByB = ZERO
        oldResNBCF = ZERO
        oldResNBCB = ZERO
        oldResNBLJF = ZERO
        oldResNBLJB = ZERO
        do i=1,NProteins
          oldResProC(i) = ZERO
          oldResProL(i) = ZERO
          oldResProCF(i) = ZERO
          oldResProLF(i) = ZERO
          oldResProCB(i) = ZERO
          oldResProLB(i) = ZERO
          oldResProDist(i) = ZERO
        enddo
        do i=1,NSolutes
          oldResSolC(i) = ZERO
          oldResSolL(i) = ZERO
          oldResSolCF(i) = ZERO
          oldResSolLF(i) = ZERO
          oldResSolCB(i) = ZERO
          oldResSolLB(i) = ZERO
          oldResSolDist(i) = ZERO
        enddo
      else if (savenrg.eq.NEW) then
        newResSvnC = ZERO
        newResSvnL = ZERO
        newResSvnCF = ZERO
        newResSvnLF = ZERO
        newResSvnCB = ZERO
        newResSvnLB = ZERO
        newResSvnDist = ZERO
	newResGCSC = ZERO
        newResGCSL = ZERO
        newResGCSCF = ZERO
        newResGCSLF = ZERO
        newResGCSCB = ZERO
        newResGCSLB = ZERO
        newResGCSDist = ZERO
        newResBnd = ZERO
        newResAng = ZERO
        newResDih = ZERO
        newResUBy = ZERO
        newResNBC = ZERO
        newResNBLJ = ZERO
        newResNrg = ZERO 
        newResNrgF = ZERO
        newResNrgB = ZERO
        newResBndF = ZERO
        newResBndB = ZERO
        newResAngF = ZERO
        newResAngB = ZERO
        newResDihF = ZERO
        newResDihB = ZERO
        newResUByF = ZERO
        newResUByB = ZERO
        newResNBCF = ZERO
        newResNBCB = ZERO
        newResNBLJF = ZERO
        newResNBLJB = ZERO
        do i=1,NProteins
          newResProC(i) = ZERO
          newResProL(i) = ZERO
          newResProCF(i) = ZERO
          newResProLF(i) = ZERO
          newResProCB(i) = ZERO
          newResProLB(i) = ZERO
          newResProDist(i) = ZERO
        enddo
        do i=1,NSolutes
          newResSolC(i) = ZERO
          newResSolL(i) = ZERO
          newResSolCF(i) = ZERO
          newResSolLF(i) = ZERO
          newResSolCB(i) = ZERO
          newResSolLB(i) = ZERO
          newResSolDist(i) = ZERO
        enddo
      endif
      
      if (iprot.le.0 .or. iprot.gt.NProteins .or. ires.le.0) return
      if (ires.gt.NResidues(iprot)) return
      
      if (DoneResNrg(iprot,ires)) then
        write(printstring,*) "Attempting to calculate the energy of ",
     .    "a residue twice? - ",iprot,ires
        call printLine(WARNING,printstring)
        return
      endif        
      
c     intra energy
      residueEnergy = residueEnergy + residueIntraEnergy(iprot,ires,
     .                bndnrg, angnrg, dihnrg, ubynrg, nbcnrg, nbljnrg,
     .                bndnrgb,angnrgb,dihnrgb,ubynrgb,nbcnrgb,nbljnrgb,
     .                bndnrgf,angnrgf,dihnrgf,ubynrgf,nbcnrgf,nbljnrgf)
     
c     inter-protein energy
      if (NProteins.gt.1) then
        do i=1,NProteins
          if (i.ne.iprot) then
            residueEnergy = residueEnergy + residueProteinEnergy(savenrg,iprot,ires,i,
     .                  icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,mindist)
            cnrg = cnrg + icnrg
            ljnrg = ljnrg + iljnrg
            cnrgf = cnrgf + icnrgf
            cnrgb = cnrgb + icnrgb
            ljnrgf = ljnrgf + iljnrgf
            ljnrgb = ljnrgb + iljnrgb
            if (savenrg.eq.OLD) then
              oldResProC(i) = icnrg
              oldResProL(i) = iljnrg
              oldResProCF(i) = icnrgf
              oldResProCB(i) = icnrgb
              oldResProLF(i) = iljnrgf
              oldResProLB(i) = iljnrgb
              oldResProDist(i) = mindist
            else if (savenrg.eq.NEW) then
              newResProC(i) = icnrg
              newResProL(i) = iljnrg
              newResProCF(i) = icnrgf
              newResProCB(i) = icnrgb
              newResProLF(i) = iljnrgf
              newResProLB(i) = iljnrgb
              newResProDist(i) = mindist 
            endif
          endif
        enddo
      endif
      
c     residue-solute energy
      if (NSolutes.gt.0) then
        do i=1,NSolutes
c         see if this is a dual topology solute, and thus whether we need to
c         scale the energies
          if (DualTopologySolute(i).gt.0) then
c           this is the low solute
            scl = getLowScale(Lambda)
            sclb = getLowScale(LambdaB)
            sclf = getLowScale(LambdaF)
          else if (DualTopologySolute(i).eq.0) then
c           this is not a dual topology solute
            scl = ONE
            sclb = ONE
            sclf = ONE
          else
c           this is the high solute
            scl = getHighScale(Lambda)
            sclb = getHighScale(LambdaB)
            sclf = getHighScale(LambdaF)
          endif         
          scl = scl*SoluteTheta(i) ! SG
          sclf = sclf*SoluteTheta(i) ! SG
          sclb = sclb*SoluteTheta(i) ! SG 
          nrg = residueSoluteEnergy(iprot,ires,i,
     .                icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,mindist,
     .                scl,sclf,sclb)
     
          if (Solutesoftcore(i)) then
            icnrg = icnrg * (scl)**CoulSoftPower
            icnrgf = icnrgf * (sclf)**CoulSoftPower
            icnrgb = icnrgb * (sclb)**CoulSoftPower
          endif         
          
          residueEnergy = residueEnergy + scl*icnrg + scl * iljnrg
          cnrg = cnrg + icnrg*scl
          ljnrg = ljnrg + iljnrg*scl
          cnrgf = cnrgf + icnrgf*sclf
          cnrgb = cnrgb + icnrgb*sclb
          ljnrgf = ljnrgf + iljnrgf*sclf
          ljnrgb = ljnrgb + iljnrgb*sclb

          if (savenrg.eq.OLD) then
            oldResSolC(i) = icnrg*scl
            oldResSolL(i) = iljnrg*scl
            oldResSolCF(i) = icnrgf*sclf
            oldResSolCB(i) = icnrgb*sclb
            oldResSolLF(i) = iljnrgf*sclf
            oldResSolLB(i) = iljnrgb*sclb
            oldResSolDist(i) = mindist
          else if (savenrg.eq.NEW) then
            newResSolC(i) = icnrg*scl
            newResSolL(i) = iljnrg*scl
            newResSolCF(i) = icnrgf*sclf
            newResSolCB(i) = icnrgb*sclb
            newResSolLF(i) = iljnrgf*sclf
            newResSolLB(i) = iljnrgb*sclb
            newResSolDist(i) = mindist
          endif
        enddo
      endif
      
c     residue-solvent energy
      if (NSolvents.gt.0) then
        residueEnergy = residueEnergy + residueSolventEnergy(savenrg,iprot,ires,
     .           icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,mindist)
        cnrg = cnrg + icnrg
        ljnrg = ljnrg + iljnrg
        cnrgf = cnrgf + icnrgf
        cnrgb = cnrgb + icnrgb
        ljnrgf = ljnrgf + iljnrgf
        ljnrgb = ljnrgb + iljnrgb
      
        if (savenrg.eq.OLD) then
          oldResSvnC = icnrg
          oldResSvnL = iljnrg
          oldResSvnCF = icnrgf
          oldResSvnCB = icnrgb
          oldResSvnLF = iljnrgf
          oldResSvnLB = iljnrgb
          oldResSvnDist = mindist
        else if (savenrg.eq.NEW) then
          newResSvnC = icnrg
          newResSvnL = iljnrg
          newResSvnCF = icnrgf
          newResSvnCB = icnrgb
          newResSvnLF = iljnrgf
          newResSvnLB = iljnrgb
          newResSvnDist = mindist
        endif
      endif
      
c     residue-GCsolute energy
      if (NGCSolutes.gt.0) then
        residueEnergy = residueEnergy + residueGCSoluteEnergy(savenrg,iprot,ires,
     .           icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,mindist)
        cnrg = cnrg + icnrg
        ljnrg = ljnrg + iljnrg
        cnrgf = cnrgf + icnrgf
        cnrgb = cnrgb + icnrgb
        ljnrgf = ljnrgf + iljnrgf
        ljnrgb = ljnrgb + iljnrgb

        if (savenrg.eq.OLD) then
          oldResGCSC = icnrg
          oldResGCSL = iljnrg
          oldResGCSCF = icnrgf
          oldResGCSCB = icnrgb
          oldResGCSLF = iljnrgf
          oldResGCSLB = iljnrgb
          oldResGCSDist = mindist
        else if (savenrg.eq.NEW) then
          newResGCSC = icnrg
          newResGCSL = iljnrg
          newResGCSCF = icnrgf
          newResGCSCB = icnrgb
          newResGCSLF = iljnrgf
          newResGCSLB = iljnrgb
          newResGCSDist = mindist
        endif
      endif
      
      nrgf = bndnrgf + angnrgf + dihnrgf + ubynrgf + nbcnrgf + nbljnrgf + cnrgf + ljnrgf
      nrgb = bndnrgb + angnrgb + dihnrgb + ubynrgb + nbcnrgb + nbljnrgb + cnrgb + ljnrgb

      if (savenrg.eq.OLD) then
        oldResBnd = bndnrg
        oldResAng = angnrg
        oldResDih = dihnrg
        oldResUBy = ubynrg
        oldResNBC = nbcnrg
        oldResNBLJ = nbljnrg
        oldResNrg = residueEnergy
        oldResNrgF = nrgf
        oldResNrgB = nrgb
        oldResBndF = bndnrgf
        oldResBndB = bndnrgb
        oldResAngF = angnrgf
        oldResAngB = angnrgb
        oldResDihF = dihnrgf
        oldResDihB = dihnrgb
        oldResUByF = ubynrgf
        oldResUByB = ubynrgb
        oldResNBCF = nbcnrgf
        oldResNBCB = nbcnrgb
        oldResNBLJF = nbljnrgf
        oldResNBLJB = nbljnrgb
      else if (savenrg.eq.NEW) then
        newResBnd = bndnrg
        newResAng = angnrg
        newResDih = dihnrg
        newResUBy = ubynrg
        newResNBC = nbcnrg
        newResNBLJ = nbljnrg
        newResNrg = residueEnergy
        newResNrgF = nrgf
        newResNrgB = nrgb
        newResBndF = bndnrgf
        newResBndB = bndnrgb
        newResAngF = angnrgf
        newResAngB = angnrgb
        newResDihF = dihnrgf
        newResDihB = dihnrgb
        newResUByF = ubynrgf
        newResUByB = ubynrgb
        newResNBCF = nbcnrgf
        newResNBCB = nbcnrgb
        newResNBLJF = nbljnrgf
        newResNBLJB = nbljnrgb
      endif    

      return
      end
