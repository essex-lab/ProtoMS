      double precision function soluteEnergy(savenrg,isol,nrgb,nrgf)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'flags.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'simparams.inc'
      include 'enums.inc'
      include 'energies.inc'
c#################################################################
c
c     This function returns the total energy of solute isol
c     with everything else (bearing in mind which things have
c     been done before).
c
c     If savenrg = 1 then this is the old energy, else if it 
c     is 2, then it is the new energy, else don't save the
c     energy!
c
c     (C) Christopher Woods
c
c     Modified by Samuel Genheden 2012-2013
c
c#################################################################

      integer i,savenrg,isol
      double precision soluteIntraEnergy,soluteSoluteEnergy
      double precision proteinSoluteEnergy,soluteSolventEnergy,soluteGCSenergy
      double precision getLowScale,getHighScale
      double precision bndnrg,angnrg,dihnrg,ubynrg,nbcnrg,nbljnrg
      double precision bndnrgf,angnrgf,dihnrgf,ubynrgf,nbcnrgf,nbljnrgf
      double precision bndnrgb,angnrgb,dihnrgb,ubynrgb,nbcnrgb,nbljnrgb
      double precision mindist,imindist,cnrg,ljnrg,icnrg,iljnrg
      double precision nrg,nrgb,nrgf
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb,icnrgf,icnrgb,iljnrgf,iljnrgb
      double precision scl,sclb,sclf


      soluteEnergy = ZERO
      mindist = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      nrgf = ZERO
      nrgb = ZERO              

      if (isol.le.0 .or. isol.gt.NSolutes) then
        return
      else if (DoneSolNrg(isol)) then
        write(printstring,*) "Already done solute energy? - ",isol
        call printLine(WARNING,printstring)
        return
      endif

c     see if this is a dual topology solute - get the scale factor if it is
      if (DualTopologySolute(isol).eq.0) then
c       not a dual topology solute
        scl = ONE
        sclf = ONE
        sclb = ONE
      else if (DualTopologySolute(isol).gt.0) then
c       this is the dual topology solute at lambda=0
        scl = getLowScale(Lambda)
        sclb = getLowScale(LambdaB)
        sclf = getLowScale(LambdaF)
      else
c       this is the dual topology solute at lambda=1
        scl = getHighScale(Lambda)
        sclb = getHighScale(LambdaB)
        sclf = getHighScale(LambdaF)
      endif  
      scl = scl*SoluteTheta(isol) ! SG
      sclf = sclf*SoluteTheta(isol) ! SG
      sclb = sclb*SoluteTheta(isol) ! SG

c     the solute intramolecular energy
      nrg = soluteIntraEnergy(isol,
     .         bndnrg, angnrg, dihnrg, ubynrg, nbcnrg, nbljnrg,
     .         bndnrgb,angnrgb,dihnrgb,ubynrgb,nbcnrgb,nbljnrgb,
     .         bndnrgf,angnrgf,dihnrgf,ubynrgf,nbcnrgf,nbljnrgf)
c     The intramolecular energy is not scaled
      soluteEnergy = soluteEnergy + nrg

      if (savenrg.eq.OLD) then
        OldSolBnd(isol) = bndnrg
        OldSolAng(isol) = angnrg
        OldSolDih(isol) = dihnrg
        OldSolUBy(isol) = ubynrg
        OldSolNBC(isol) = nbcnrg
        OldSolNBLJ(isol) = nbljnrg
        OldSolBndF(isol) = bndnrgf
        OldSolAngF(isol) = angnrgf
        OldSolDihF(isol) = dihnrgf
        OldSolUByF(isol) = ubynrgf
        OldSolNBCF(isol) = nbcnrgf
        OldSolNBLJF(isol) = nbljnrgf
        OldSolBndB(isol) = bndnrgb
        OldSolAngB(isol) = angnrgb
        OldSolDihB(isol) = dihnrgb
        OldSolUByB(isol) = ubynrgb
        OldSolNBCB(isol) = nbcnrgb
        OldSolNBLJB(isol) = nbljnrgb
      else if (savenrg.eq.NEW) then
        NewSolBnd(isol) = bndnrg
        NewSolAng(isol) = angnrg
        NewSolDih(isol) = dihnrg
        NewSolUBy(isol) = ubynrg
        NewSolNBC(isol) = nbcnrg
        NewSolNBLJ(isol) = nbljnrg
        NewSolBndF(isol) = bndnrgf
        NewSolAngF(isol) = angnrgf
        NewSolDihF(isol) = dihnrgf
        NewSolUByF(isol) = ubynrgf
        NewSolNBCF(isol) = nbcnrgf
        NewSolNBLJF(isol) = nbljnrgf
        NewSolBndB(isol) = bndnrgb
        NewSolAngB(isol) = angnrgb
        NewSolDihB(isol) = dihnrgb
        NewSolUByB(isol) = ubynrgb
        NewSolNBCB(isol) = nbcnrgb
        NewSolNBLJB(isol) = nbljnrgb
      endif        

c     inter-solute energy
      if (NSolutes.gt.1) then
        do i=1,NSolutes
          if (i.ne.isol) then
c           the scaling of solute-solute energies to account for 
c           multi-copy solutes is performed in soluteSoluteEnergy 
c           due to the complexity of potentially having
c           multiple solutes involved in dual topology sampling
            soluteEnergy = soluteEnergy + soluteSoluteEnergy(isol,i,
     .             icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
            cnrg = cnrg + icnrg
            ljnrg = ljnrg + iljnrg
            cnrgf = cnrgf + icnrgf
            cnrgb = cnrgb + icnrgb
            ljnrgf = ljnrgf + iljnrgf
            ljnrgb = ljnrgb + iljnrgb
            if (savenrg.eq.OLD) then
              OldSolSolC(isol,i) = icnrg
              OldSolSolLJ(isol,i) = iljnrg
              OldSolSolC(i,isol) = icnrg
              OldSolSolLJ(i,isol) = iljnrg
              OldSolSolCF(isol,i) = icnrgf
              OldSolSolLJF(isol,i) = iljnrgf
              OldSolSolCF(i,isol) = icnrgf
              OldSolSolLJF(i,isol) = iljnrgf
              OldSolSolCB(isol,i) = icnrgb
              OldSolSolLJB(isol,i) = iljnrgb
              OldSolSolCB(i,isol) = icnrgb
              OldSolSolLJB(i,isol) = iljnrgb
              OldSolSolDist(isol,i) = imindist
              OldSolSolDist(i,isol) = imindist
            else if (savenrg.eq.NEW) then
              NewSolSolC(isol,i) = icnrg
              NewSolSolLJ(isol,i) = iljnrg
              NewSolSolC(i,isol) = icnrg
              NewSolSolLJ(i,isol) = iljnrg
              NewSolSolCF(isol,i) = icnrgf
              NewSolSolLJF(isol,i) = iljnrgf
              NewSolSolCF(i,isol) = icnrgf
              NewSolSolLJF(i,isol) = iljnrgf
              NewSolSolCB(isol,i) = icnrgb
              NewSolSolLJB(isol,i) = iljnrgb
              NewSolSolCB(i,isol) = icnrgb
              NewSolSolLJB(i,isol) = ljnrgb
              NewSolSolDist(isol,i) = imindist
              NewSolSolDist(i,isol) = imindist
            endif      
          endif
        enddo
      endif
      
c     solute-protein energy
      if (NProteins.gt.0) then
        do i=1,NProteins
          nrg = proteinSoluteEnergy(SAVENRG,i,isol,
     .              icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist,
     .              scl,sclf,sclb)
c$$$c         If the soft-core is on, the coulombic is scaled further here
c$$$          if (Solutesoftcore(isol)) then
c$$$            icnrg = icnrg * (scl)**CoulSoftPower
c$$$            icnrgf = icnrgf * (sclf)**CoulSoftPower
c$$$            icnrgb = icnrgb * (sclb)**CoulSoftPower
c$$$          endif                
          cnrg = cnrg + icnrg*scl
          ljnrg = ljnrg + iljnrg*scl
          cnrgf = cnrgf + icnrgf*sclf
          cnrgb = cnrgb + icnrgb*sclb
          ljnrgf = ljnrgf + iljnrgf*sclf
          ljnrgb = ljnrgb + iljnrgb*sclb
!         SG
          soluteEnergy = soluteEnergy + scl*icnrg + scl*iljnrg

          if (savenrg.eq.OLD) then
            OldProSolC(i,isol) = icnrg*scl
            OldProSolLJ(i,isol) = iljnrg*scl
            OldProSolCF(i,isol) = icnrgf*sclf
            OldProSolLJF(i,isol) = iljnrgf*sclf
            OldProSolCB(i,isol) = icnrgb*sclb
            OldProSolLJB(i,isol) = iljnrgb*sclb
            OldProSolDist(i,isol) = imindist
          else if (savenrg.eq.NEW) then
            NewProSolC(i,isol) = icnrg*scl
            NewProSolLJ(i,isol) = iljnrg*scl
            NewProSolCF(i,isol) = icnrgf*sclf
            NewProSolLJF(i,isol) = iljnrgf*sclf
            NewProSolCB(i,isol) = icnrgb*sclb
            NewProSolLJB(i,isol) = iljnrgb*sclb
            NewProSolDist(i,isol) = imindist
          endif      
        enddo
      endif
      
c    solute-gcs energy
      if (NGCSolutes.gt.0) then
       nrg = soluteGCSEnergy(savenrg,isol,
     .               icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist,
     .               scl,sclf,sclb)

!       SG
c       If the soft-core is on, the coulombic is scaled further here
        if (Solutesoftcore(isol)) then
          icnrg = icnrg * (scl)**CoulSoftPower
          icnrgf = icnrgf * (sclf)**CoulSoftPower
          icnrgb = icnrgb * (sclb)**CoulSoftPower
        endif   
      
        cnrg = cnrg + icnrg*scl
        ljnrg = ljnrg + iljnrg*scl
        cnrgf = cnrgf + icnrgf*sclf
        cnrgb = cnrgb + icnrgb*sclb
        ljnrgf = ljnrgf + iljnrgf*sclf
        ljnrgb = ljnrgb + iljnrgb*sclb
        soluteEnergy = soluteEnergy + icnrg*scl + iljnrg*scl

        if (savenrg.eq.OLD) then
          OldSolGCSC(isol) = icnrg*scl
          OldSolGCSLJ(isol) = iljnrg*scl
          OldSolGCSCF(isol) = icnrgf*sclf
          OldSolGCSLJF(isol) = iljnrgf*sclf
          OldSolGCSCB(isol) = icnrgb*sclb
          OldSolGCSLJB(isol) = iljnrgb*sclb
        else if (savenrg.eq.NEW) then
          NewSolGCSC(isol) = icnrg*scl
          NewSolGCSLJ(isol) = iljnrg*scl
          NewSolGCSCF(isol) = icnrgf*sclf
          NewSolGCSLJF(isol) = iljnrgf*sclf
          NewSolGCSCB(isol) = icnrgb*sclb
          NewSolGCSLJB(isol) = iljnrgb*sclb
        endif
      endif

      
c     solute-solvent energy
      if (NSolvents.gt.0) then
        nrg = soluteSolventEnergy(savenrg,isol,
     .               icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist,
     .               scl,sclf,sclb)
c       If the soft-core is on, the coulombic is scaled further here
c$$$        if (Solutesoftcore(isol)) then
c$$$          icnrg = icnrg * (scl)**CoulSoftPower
c$$$          icnrgf = icnrgf * (sclf)**CoulSoftPower
c$$$          icnrgb = icnrgb * (sclb)**CoulSoftPower
c$$$        endif         
        cnrg = cnrg + icnrg*scl
        ljnrg = ljnrg + iljnrg*scl
        cnrgf = cnrgf + icnrgf*sclf
        cnrgb = cnrgb + icnrgb*sclb
        ljnrgf = ljnrgf + iljnrgf*sclf
        ljnrgb = ljnrgb + iljnrgb*sclb
        soluteEnergy = soluteEnergy + scl*icnrg + scl*iljnrg
          
        if (savenrg.eq.OLD) then
          OldSolSvnC(isol) = icnrg*scl
          OldSolSvnLJ(isol) = iljnrg*scl
          OldSolSvnCF(isol) = icnrgf*sclf
          OldSolSvnLJF(isol) = iljnrgf*sclf
          OldSolSvnCB(isol) = icnrgb*sclb
          OldSolSvnLJB(isol) = iljnrgb*sclb
          OldSolSvnDist(isol) = imindist
c         if this is the solute used for preferential sampling,
c         then save the preferential sampling parameters
          if (isol.eq.SvnPrefSol .and. PrefSample) then
            do i=1,NSolvents
              SvnProb(i) = SvnNewProb(i)
            enddo
            SvnSumProb = SvnNewSumProb
            SvnMaxProb = SvnNewMaxProb
          endif
        else if (savenrg.eq.NEW) then
          NewSolSvnC(isol) = icnrg*scl
          NewSolSvnLJ(isol) = iljnrg*scl
          NewSolSvnCF(isol) = icnrgf*sclf
          NewSolSvnLJF(isol) = iljnrgf*sclf
          NewSolSvnCB(isol) = icnrgb*sclb
          NewSolSvnLJB(isol) = iljnrgb*sclb
          NewSolSvnDist(isol) = imindist
        endif      
      endif

      nrgf = (bndnrgf+angnrgf+dihnrgf+ubynrgf+nbcnrgf+nbljnrgf) +cnrgf+ljnrgf
      nrgb = (bndnrgb+angnrgb+dihnrgb+ubynrgb+nbcnrgb+nbljnrgb) +cnrgb+ljnrgb

      if (savenrg.eq.OLD) then
        oldSoluteEnergy(isol) = soluteEnergy
        oldSoluteEnergyF(isol) = nrgf       
        oldSoluteEnergyB(isol) = nrgb
      else if (savenrg.eq.NEW) then
        newSoluteEnergy(isol) = soluteEnergy
        newSoluteEnergyF(isol) = nrgf
        newSoluteEnergyB(isol) = nrgb
      endif
c      write(printstring,*) 'solNrg ',soluteEnergy,' F ',nrgf,' B ',nrgb
c      call printLine(DEBUG,printstring)


      return
      end
