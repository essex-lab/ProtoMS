      
      double precision function soluteGCSEnergy(savecut,isol,
     .            cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist,
     .            scl,sclf,sclb)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'flags.inc'
c################################################################
c
c     This function calculates the total interaction energy
c     between solute isol, and all the gc-solutes
c
c     (C) MSB   April 2012
c
c################################################################
      include 'printer.inc'
      integer isol,savecut
      integer i
      double precision cnrg,ljnrg,nrg,icnrg,iljnrg,mindist,imindist,inrg
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb,icnrgf,icnrgb,iljnrgf,iljnrgb
      double precision soluteSingleGCSoluteEnergy,idel
      double precision scl,sclf,sclb
      soluteGCSEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      nrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      mindist = 999999.0

      if (isol.gt.NSolutes  .or. NSolutes.le.0 .or. NGCSolutes.le.0) return

      if (IgnoreSolute(isol)) return   
    
c     return if this solute has been done before
      if (DoneSolNrg(isol)) return
      do i=1,NGCSolutes
        inrg = soluteSingleGCSoluteEnergy(isol,i,
     .         icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist,
     .         scl,sclf,sclb)
        nrg = nrg + inrg
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
      if (isol.eq.SvnPrefSol .and. PrefSample) then
c       initialise the sum and max sampling
        SvnNewSumProb = ZERO
        SvnNewMaxProb = ZERO

c       calculate the new sum of sampling
        do i=1,NSolvents
          SvnNewProb(i) = ONE / (SvnNewProb(i) + SvnProbParam)
          SvnNewSumProb = SvnNewSumProb + SvnNewProb(i)
        enddo

c       normalise the SvnNewProb... - while doing so, watch 
c       for the maximum probability..
        idel = ONE / SvnNewSumProb
        do i=1,NSolvents
          SvnNewProb(i) = idel * SvnNewProb(i)
          if (SvnNewProb(i).gt.SvnNewMaxProb) SvnNewMaxProb = SvnNewProb(i)
        enddo
      endif

      soluteGCSEnergy = nrg

      return
      end
