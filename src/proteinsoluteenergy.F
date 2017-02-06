      double precision function proteinSoluteEnergy(savecut,iprot,isol,cnrg,ljnrg,
     .  cnrgb,ljnrgb,cnrgf,ljnrgf,mindist,scl,sclf,sclb)
      include 'dimensions.inc'
      include 'proteins.inc'
      include 'simparams.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'enums.inc'
c#############################################################
c
c     This function returns the total interaction energy between
c     protein iprot and solute isol
c
c     (C) Christopher Woods      29 January 2002
c
c#############################################################

      integer iprot,isol,savecut
      double precision cnrg,ljnrg,nrg,icnrg,iljnrg
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb,icnrgf,icnrgb,iljnrgf,iljnrgb
      double precision mindist,imindist,sclfac
      integer i,j,nres
      double precision residueSoluteEnergy
      double precision scl,sclf,sclb
      double precision idel
      proteinSoluteEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      nrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      mindist = 99000.0
      
      if (NProteins.le.0 .or. iprot.gt.NProteins .or. NSolutes.le.0
     .         .or. isol.gt.NSolutes) return

c      PRINT*, "I AM",isol 
      if (IgnoreSolute(isol)) return
          
      nres = NResidues(iprot)
      if (nres.le.0) return
      
      do i=1,nres
        nrg = nrg + residueSoluteEnergy(iprot,i,isol,
     .               icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist,
     .               scl,sclf,sclb)
        cnrg = cnrg + icnrg
        ljnrg = ljnrg + iljnrg
        cnrgf = cnrgf + icnrgf
        cnrgb = cnrgb + icnrgb
        ljnrgf = ljnrgf + iljnrgf
        ljnrgb = ljnrgb + iljnrgb
        if (imindist.lt.mindist) mindist = imindist
      enddo

c     Sort out the preferential sampling - we need to normalise
c     all of the probabilities, and find the maximum probability
      if (isol.eq.ProPrefSol .and. ProPrefSample .and. .not. doDerivative) then
c       initialise the sum and max sampling
        ProNewSumProb = ZERO
        ProNewMaxProb = ZERO

c       calculate the new sum of sampling
        do i=1,NProteins
           if (i.eq.iprot) then
c     only take the new probability values for the protein we've calculated for
              do j=1,NResidues(iprot)
                 ProNewProb(i,j) = ONE / (ProNewProb(i,j) + ProProbParam)
                 ProNewSumProb = ProNewSumProb + ProNewProb(i,j)
              enddo
           else
c     take the unnormalised values for the other proteins
              do j=1,NResidues(i)
                 ProNewProb(i,j) = ProProb(i,j) * ProSumProb
                 ProNewSumProb = ProNewSumProb + ProNewProb(i,j)
              enddo
           endif
        enddo

c       normalise the ProNewProb... - while doing so, watch 
c       for the maximum probability..
        idel = ONE / ProNewSumProb
        do i=1,NProteins
           do j=1,NResidues(i)
              ProNewProb(i,j) = idel * ProNewProb(i,j)
              if (ProNewProb(i,j).gt.ProNewMaxProb) ProNewMaxProb = ProNewProb(i,j)
           enddo
        enddo
      endif
      
      
      if (cuttype.eq.MOLECULEMOLECULE) then
c       performing protein-molecule solute-molecule cutting
        if (mindist.lt.cut2) then
c         within the cutoff distance
          if (mindist.gt.feather2) then
c           however, in the feather region
            sclfac = featherfac * (cut2-mindist)
            nrg = nrg * sclfac
            cnrg = cnrg * sclfac
            ljnrg = ljnrg * sclfac
            cnrgf = cnrgf * sclfac
            ljnrgf = ljnrgf * sclfac
            cnrgb = cnrgb * sclfac
            ljnrgb = ljnrgb * sclfac
          endif
        else
          nrg = ZERO
          cnrg = ZERO
          ljnrg = ZERO
          cnrgf = ZERO
          cnrgb = ZERO
          ljnrgf = ZERO
          ljnrgb = ZERO
        endif
      endif
      
      proteinSoluteEnergy = nrg
      
      return
      end
