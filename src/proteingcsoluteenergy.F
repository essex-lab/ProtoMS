      double precision function proteinGCSoluteEnergy(iprot,
     .                   cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'enums.inc'
      include 'constants.inc'
      include 'proteins.inc'
c#####################################################################
c
c     This function returns the protein-GCS energy between
c     protein iprot, and all the GCS molecules
c     It also returns the coulomb and LJ components in cnrg and ljnrg
c
c     (C) Christopher Woods    29 January 2002
c       MSB Sep 2012
c#####################################################################

      integer iprot
      double precision nrg,cnrg,ljnrg,icnrg,iljnrg
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb,icnrgf,icnrgb,iljnrgf,iljnrgb
      integer i,nres
      double precision residueGCSoluteEnergy
      double precision mindist,imindist

      proteinGCSoluteEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      mindist = 99000.0
      if (NGCSolutes.le.0 .or. NProteins.le.0 .or. iprot.gt.NProteins) return
      
      nrg = ZERO
      nres = NResidues(iprot)
      
      do i=1,nres
        nrg = nrg + residueGCSoluteEnergy(DONTCARE,iprot,i,
     .              icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
        cnrg = cnrg + icnrg
        ljnrg = ljnrg + iljnrg
        cnrgf = cnrgf + icnrgf
        cnrgb = cnrgb + icnrgb
        ljnrgf = ljnrgf + iljnrgf
        ljnrgb = ljnrgb + iljnrgb
        if (imindist.lt.mindist) mindist = imindist
      enddo
      
      proteinGCSoluteEnergy = nrg

      return
      
      end
