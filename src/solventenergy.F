      double precision function solventEnergy(savenrg,isvn,nrgb,nrgf)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'flags.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'solvents.inc'
      include 'simparams.inc'
      include 'enums.inc'
      include 'energies.inc'
c ### delete later below
      include 'parallel.inc'
c##############################################################
c
c     This function returns the total energy of solvent isvn,
c     (bearing in mind what has come before)
c
c     If savenrg = 1 then this is the old energy, else if it 
c     is 2, then it is the new energy, else don't save the
c     energy!
c
c     (C) Christopher Woods  14 February 2002
c
c#################################################################

      integer i,j,savenrg,isvn,nres
      double precision cnrg,ljnrg,mindist,icnrg,iljnrg,imindist,inrg
      double precision iicnrg,iiljnrg,iicnrgf,iicnrgb,iiljnrgf,iiljnrgb
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb,icnrgf,icnrgb,iljnrgf,iljnrgb
      double precision nrg,scl,sclf,sclb
      double precision getHighScale,getLowScale
      double precision solventSolventEnergy,soluteSingleSolventEnergy
      double precision residueSingleSolventEnergy,GCsoluteSingleSolventEnergy
      double precision nrgb,nrgf,idel

      solventEnergy = ZERO
      nrgf = ZERO
      nrgb = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      ljnrgf = ZERO
      cnrgb = ZERO
      ljnrgb = ZERO
      
      if (NSolvents.le.0 .or. isvn.gt.MAXSOLVENTS) return
      
      if (DoneSvnNrg(isvn)) then
        write(printstring,*) "Have we already done this solvent? - ",isvn
        call printLine(WARNING,printstring)
        return
      endif
      
      inrg = solventSolventEnergy(savenrg,isvn,
     .               icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,mindist)
      solventEnergy = solventEnergy + inrg
      cnrg = icnrg
      ljnrg = iljnrg
      cnrgf = icnrgf
      cnrgb = icnrgb
      ljnrgf = iljnrgf
      ljnrgb = iljnrgb
      
      if (savenrg.eq.NEW) then
        NewIndSvnSvnC = icnrg
        NewIndSvnSvnLJ = iljnrg
        NewIndSvnSvnCF = icnrgf
        NewIndSvnSvnLJF = iljnrgf
        NewIndSvnSvnCB = icnrgb
        NewIndSvnSvnLJB = iljnrgb
      else if (savenrg.eq.OLD) then
        OldIndSvnSvnC = icnrg
        OldIndSvnSvnLJ = iljnrg
        OldIndSvnSvnCF = icnrgf
        OldIndSvnSvnLJF = iljnrgf
        OldIndSvnSvnCB = icnrgb
        OldIndSvnSvnLJB = iljnrgb
      endif
      
      if (NProteins.gt.0) then
        do i=1,NProteins
          nres = NResidues(i)
          iicnrg = ZERO
          iiljnrg = ZERO
          iicnrgf = ZERO
          iiljnrgf = ZERO
          iicnrgb = ZERO
          iiljnrgb = ZERO
          do j=1,nres
            solventEnergy = solventEnergy + residueSingleSolventEnergy(i,j,isvn,
     .            icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
            cnrg = cnrg + icnrg
            ljnrg = ljnrg + iljnrg
            cnrgf = cnrgf + icnrgf
            cnrgb = cnrgb + icnrgb
            ljnrgf = ljnrgf + iljnrgf
            ljnrgb = ljnrgb + iljnrgb
            iicnrg = iicnrg + icnrg
            iiljnrg = iiljnrg + iljnrg
            iicnrgf = iicnrgf + icnrgf
            iicnrgb = iicnrgb + icnrgb
            iiljnrgf = iiljnrgf + iljnrgf
            iiljnrgb = iiljnrgb + iljnrgb
            if (imindist.lt.mindist) mindist = imindist
          enddo
          if (savenrg.eq.NEW) then
            NewIndSvnProC(i) = iicnrg
            NewIndSvnProLJ(i) = iiljnrg
            NewIndSvnProCF(i) = iicnrgf
            NewIndSvnProLJF(i) = iiljnrgf
            NewIndSvnProCB(i) = iicnrgb
            NewIndSvnProLJB(i) = iiljnrgb
          else if (savenrg.eq.OLD) then
            OldIndSvnProC(i) = iicnrg
            OldIndSvnProLJ(i) = iiljnrg
            OldIndSvnProCF(i) = iicnrgf
            OldIndSvnProLJF(i) = iiljnrgf
            OldIndSvnProCB(i) = iicnrgb
            OldIndSvnProLJB(i) = iiljnrgb
          endif
        enddo
      endif
      
      if (NGCSolutes.gt.0) then
        do i=1,NGCSolutes
        solventEnergy = solventEnergy + GCsoluteSingleSolventEnergy(i,isvn,
     .                    icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
            cnrg = cnrg + icnrg
            ljnrg = ljnrg + iljnrg
            cnrgf = cnrgf + icnrgf
            cnrgb = cnrgb + icnrgb
            ljnrgf = ljnrgf + iljnrgf
            ljnrgb = ljnrgb + iljnrgb
          if (savenrg.eq.NEW) then
            NewIndSvnGCSC(i) = icnrg
            NewIndSvnGCSLJ(i) = iljnrg
            NewIndSvnGCSCF(i) = icnrgf
            NewIndSvnGCSLJF(i) = iljnrgf
            NewIndSvnGCSCB(i) = icnrgb
            NewIndSvnGCSLJB(i) = iljnrgb
          else if (savenrg.eq.OLD) then
            OldIndSvnGCSC(i) = icnrg
            OldIndSvnGCSLJ(i) = iljnrg
            OldIndSvnGCSCF(i) = icnrgf
            OldIndSvnGCSLJF(i) = iljnrgf
            OldIndSvnGCSCB(i) = icnrgb
            OldIndSvnGCSLJB(i) = iljnrgb
          endif
        enddo
      endif
      
      if (NSolutes.gt.0) then
        do i=1,NSolutes
c         see whether or not we need to scale the energy
          if (DualTopologySolute(i).gt.0) then
c           this is the low dual topology solute
            scl = getLowScale(Lambda)
            sclb = getLowScale(LambdaB)
            sclf = getLowScale(LambdaF)
          else if (DualTopologySolute(i).eq.0) then
c           this is not a dual topology solute
            scl = ONE
            sclb = ONE
            sclf = ONE
          else
c           this is the high dual topology solute
            scl = getHighScale(Lambda)
            sclb = getHighScale(LambdaB)
            sclf = getHighScale(LambdaF)
          endif
          scl = scl*SoluteTheta(i) ! SG
          sclf = sclf*SoluteTheta(i) ! SG
          sclb = sclb*SoluteTheta(i) ! SG 
          nrg = soluteSingleSolventEnergy(i,isvn,
     .        icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist,
     .        scl,sclf,sclb)

          if (Solutesoftcore(i)) then
            icnrg = icnrg * (scl)**CoulSoftPower
            icnrgf = icnrgf * (sclf)**CoulSoftPower
            icnrgb = icnrgb * (sclb)**CoulSoftPower
          endif       
          solventEnergy = solventEnergy + scl*icnrg + scl*iljnrg

          cnrg = cnrg + icnrg*scl
          ljnrg = ljnrg + iljnrg*scl
          cnrgf = cnrgf + icnrgf*sclf
          cnrgb = cnrgb + icnrgb*sclb
          ljnrgf = ljnrgf + iljnrgf*sclf
          ljnrgb = ljnrgb + iljnrgb*sclb
          if (savenrg.eq.NEW) then
            NewIndSvnSolC(i) = icnrg*scl
            NewIndSvnSolLJ(i) = iljnrg*scl
            NewIndSvnSolCF(i) = icnrgf*sclf
            NewIndSvnSolLJF(i) = iljnrgf*sclf
            NewIndSvnSolCB(i) = icnrgb*sclb
            NewIndSvnSolLJB(i) = iljnrgb*sclb
c           if this is a new conformation, and we are interacting with
c           the preferential solute, then update the solvent probabilities
            if (i.eq.SvnPrefSol .and. PrefSample) then
              SvnNewSumProb = ZERO
              SvnNewMaxProb = ZERO
c             unnormalise the probabilities! - sum together as well
              do j=1,NSolvents
                if (j.ne.isvn) then
                  SvnNewProb(j) = SvnProb(j) * SvnSumProb
                  SvnNewSumProb = SvnNewSumProb + SvnNewProb(j)
                endif                  
              enddo
c             calculate the probability of this solvent
              SvnNewProb(isvn) = ONE / (SvnNewProb(isvn) + SvnProbParam)
c             add this back onto the sum, and normalise the probabilities
c             - also watch for the maximum probability
              SvnNewSumProb = SvnNewSumProb + SvnNewProb(isvn)
              idel = ONE / SvnNewSumProb
              do j=1,NSolvents              
                SvnNewProb(j) = idel * SvnNewProb(j)
                if (SvnNewProb(j).gt.SvnNewMaxProb) SvnNewMaxProb = SvnNewProb(j)
              enddo
            endif
          else if (savenrg.eq.OLD) then
            OldIndSvnSolC(i) = icnrg*scl
            OldIndSvnSolLJ(i) = iljnrg*scl
            OldIndSvnSolCF(i) = icnrgf*sclf
            OldIndSvnSolLJF(i) = iljnrgf*sclf
            OldIndSvnSolCB(i) = icnrgb*sclb
            OldIndSvnSolLJB(i) = iljnrgb*sclb
          endif
        enddo
      endif
      
      if (savenrg.eq.NEW) then
        NewIndSvnC = cnrg
        NewIndSvnLJ = ljnrg        
        NewIndSvnCF = cnrgf
        NewIndSvnLJF = ljnrgf
        NewIndSvnCB = cnrgb
        NewIndSvnLJB = ljnrgb
      else if (savenrg.eq.OLD) then
        OldIndSvnC = cnrg
        OldIndSvnLJ = ljnrg  
        OldIndSvnCF = cnrgf
        OldIndSvnLJF = ljnrgf
        OldIndSvnCB = cnrgb
        OldIndSvnLJB = ljnrgb
      endif
      
      nrgf = cnrgf + ljnrgf
      nrgb = cnrgb + ljnrgb

      if (myid.eq.1) write(6,*) "For solventenergy NewIndSvnSolC(1)",NewIndSvnSolC(1)
      
      return
      end
