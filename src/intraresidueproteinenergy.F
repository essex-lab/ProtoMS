      double precision function intraResidueProteinEnergy(iprot,ires,
     .            cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'proteins.inc'
c####################################################################
c
c     This function calculates the non-bonded interaction energy between
c     residue ires with the rest of protein iprot.
c     It also returns the coulomb and lj components in cnrg and ljnrg...
c
c     (C) Christopher Woods  25 January 2002
c
c####################################################################

      integer iprot,ires
      double precision cnrg,ljnrg,icnrg,iljnrg,nrg,inrg
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb,icnrgf,icnrgb,iljnrgf,iljnrgb
      integer i,inres
      double precision residueResidueEnergy
      double precision mindist,imindist

      intraResidueProteinEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      mindist = 99000.0

      if (ires.le. 0 .or. iprot.le.0 .or. iprot.gt.MAXPROTEINS) return
      inres = NResidues(iprot)
      if (ires.gt.inres .or. ires.le.0) return
      
      if (DoneResNrg(iprot,ires)) return
      
      nrg = ZERO
      
c     loop over all of the residues in the protein, except this one!
      do i=1,inres      
        if (i.ne.ires) then
          inrg = residueResidueEnergy(iprot,ires,i,
     .                icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)

          nrg = nrg + inrg
          cnrg = cnrg + icnrg
          ljnrg = ljnrg + iljnrg
          cnrgf = cnrgf + icnrgf
          cnrgb = cnrgb + icnrgb
          ljnrgf = ljnrgf + iljnrgf
          ljnrgb = ljnrgb + iljnrgb
          
          if (imindist.lt.mindist) mindist = imindist
        endif
      enddo

      intraResidueProteinEnergy = nrg

      return
      
      end
