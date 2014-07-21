      double precision function proteinProteinEnergy(iprot1,iprot2,cnrg,ljnrg,
     .                          cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'constants.inc'
      include 'proteins.inc'
c#############################################################
c
c     This function returns the inter-protein
c     energy (sum of interprotein LJ and electostatic
c     terms) calculated from scratch...
c     The coulomb and LJ components are retuned in cnrg and ljnrg...
c
c     (C) Christopher Woods   23 January 2002
c
c#############################################################
      
      integer iprot1,iprot2,nres1,nres2
      integer i
      double precision cnrg,ljnrg,icnrg,iljnrg,nrg,cnrgb,cnrgf,ljnrgb,ljnrgf
      double precision icnrgf,icnrgb,iljnrgf,iljnrgb
      double precision mindist,imindist
      double precision residueProteinEnergy

      proteinProteinEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgb = ZERO
      cnrgf = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      nrg = ZERO
      mindist = 99000.0

      if (iprot1.eq.iprot2) return
      if (iprot1.le.0 .or. iprot1.gt.NProteins
     .      .or. iprot2.le.0  .or. iprot2.gt.NProteins)  return
      
      nres1 = NResidues(iprot1)
      nres2 = NResidues(iprot2)
      
      if (nres1.le.0  .or. nres2.le.0) return

      do i=1,nres1
        if (.not.DoneResNrg(iprot1,i)) then
c         this would never be used to perform moves - thus we don't
c         care whether this is old or new
          nrg = nrg + residueProteinEnergy(DONTCARE,iprot1,i,iprot2,
     .                  icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
          cnrg = cnrg + icnrg
          ljnrg = ljnrg + iljnrg
          cnrgf = cnrgf + icnrgf
          cnrgb = cnrgb + icnrgb
          ljnrgf = ljnrgf + iljnrgf
          ljnrgb = ljnrgb + iljnrgb
          if (imindist.lt.mindist) mindist = imindist
        endif
      enddo
      
      proteinProteinEnergy = nrg
      
      return 
      
      end
