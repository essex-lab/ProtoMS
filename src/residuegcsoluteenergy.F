      double precision function residueGCSoluteEnergy(savecut,iprot,ires,
     .             cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'constants.inc'
      include 'proteins.inc'
c#########################################################################
c
c     This calculates the interaction energy between residue ires of
c     protein iprot, and all of the GC-solutes...
c
c     (C) Michael Bodnarchuk     08 December 2009
c
c#########################################################################

      integer iprot,ires,savecut
      double precision cnrg,ljnrg,nrg,icnrg,iljnrg
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb,icnrgf,icnrgb,iljnrgf,iljnrgb
      integer i
      double precision mindist,imindist
      double precision residueSingleGCsoluteEnergy

      residueGCSoluteEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      nrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      mindist = 999000.0

      if (NProteins.le.0 .or. iprot.gt.NProteins .or. ires.le.0
     .        .or. NGCSolutes.le.0) return
     
      if (ires .gt. NResidues(iprot)) return

c     return if the residue has been done before
      if (DoneResNrg(iprot,ires)) return

      do i=1,NGCSolutes
        nrg = nrg + residueSingleGCsoluteEnergy(iprot,ires,i,
     .         icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
        cnrg = cnrg + icnrg
        ljnrg = ljnrg + iljnrg
        cnrgf = cnrgf + icnrgf
        cnrgb = cnrgb + icnrgb
        ljnrgf = ljnrgf + iljnrgf
        ljnrgb = ljnrgb + iljnrgb
        if (imindist.lt.mindist) mindist = imindist
      enddo
      
      residueGCSoluteEnergy = nrg
      return
      end
