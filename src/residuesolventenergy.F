      double precision function residueSolventEnergy(savecut,iprot,ires,
     .             cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'constants.inc'
      include 'proteins.inc'
c#########################################################################
c
c     This calculates the interaction energy between residue ires of
c     protein iprot, and all of the solvents...
c
c     (C) Christopher Woods     29 January 2002
c
c#########################################################################

      integer iprot,ires,savecut
      double precision cnrg,ljnrg,nrg,icnrg,iljnrg
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb,icnrgf,icnrgb,iljnrgf,iljnrgb
      integer i
      double precision mindist,imindist
      double precision residueSingleSolventEnergy

      residueSolventEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      nrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      mindist = 999000.0

      if (NProteins.le.0 .or. iprot.gt.NProteins .or. ires.le.0
     .        .or. NSolvents.le.0) return
     
      if (ires .gt. NResidues(iprot)) return

c     return if the residue has been done before
      if (DoneResNrg(iprot,ires)) return

c!$OMP PARALLEL DO DEFAULT(NONE)
c!$OMP&SHARED(iprot,ires)
c!$OMP&PRIVATE(icnrg,iljnrg,icnrgf,icnrgb,iljnrgf,iljnrgb,imindist)
c!$OMP&REDUCTION(+:nrg,cnrg,ljnrg,cnrgf,cnrgb,ljnrgf,ljnrgb)
c!$OMP&REDUCTION(min:mindist)
      do i=1,NSolvents
        nrg = nrg + residueSingleSolventEnergy(iprot,ires,i,
     .         icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
        cnrg = cnrg + icnrg
        ljnrg = ljnrg + iljnrg
        cnrgf = cnrgf + icnrgf
        cnrgb = cnrgb + icnrgb
        ljnrgf = ljnrgf + iljnrgf
        ljnrgb = ljnrgb + iljnrgb
        if (imindist.lt.mindist) mindist = imindist
      enddo
c!$OMP END PARALLEL DO
      
      residueSolventEnergy = nrg
      return
      end
