      double precision function residueProteinEnergy(savecut,iprot1,ires,iprot2,
     .             cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
      include 'dimensions.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'simparams.inc'
      include 'flags.inc'
      include 'enums.inc'
c####################################################################
c
c     This function calculates the non-bonded interaction energy between
c     residue ires of protein iprot1, and the whole of iprot2.
c     It also returns the coulomb and lj components in cnrg and ljnrg...
c
c     (C) Christopher Woods  25 January 2002
c
c####################################################################

      integer iprot1,ires,iprot2,i,inres,savecut
      double precision cnrg,ljnrg,icnrg,iljnrg,nrg
      double precision cnrgf,ljnrgf,icnrgf,iljnrgf
      double precision cnrgb,ljnrgb,icnrgb,iljnrgb
      double precision mindist,imindist,sclfac
      double precision intraResidueProteinEnergy
      double precision interResidueResidueEnergy

c#################################################
C     CVS keyword - used to tag the binary
c      character*80 id,name,date
c      name = "$Author: chris $"
c      date = "$Date: 2002/08/12 20:11:06 $"
c      id = "$Id: residueproteinenergy.F,v 1.12 2002/08/12 20:11:06 chris Exp $"
c#################################################

      residueProteinEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      nrg = ZERO
      mindist = 99999.0

      if (iprot1.eq.iprot2) then
c       use a different function if calculating the intramolecular
c       residue-protein energy (since need non-bonded/bonded interactions...)
        residueProteinEnergy = intraResidueProteinEnergy(iprot1,ires,
     .               cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
        return
      endif

      if (DoneResNrg(iprot1,ires)) return

      if (iprot1.le.0 .or. iprot1.gt.MAXPROTEINS
     .      .or. iprot2.le.0  .or. iprot2.gt.MAXPROTEINS)  return
      
      inres = NResidues(iprot2)
      if (ires.gt.NResidues(iprot1) .or. inres.le.0) return
      
      do i=1,inres
        nrg = nrg + interResidueResidueEnergy(iprot1,ires,iprot2,i,
     .          icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
        cnrg = cnrg + icnrg
        ljnrg = ljnrg + iljnrg
        cnrgf = cnrgf + icnrgf
        cnrgb = cnrgb + icnrgb
        ljnrgf = ljnrgf + iljnrgf
        ljnrgb = ljnrgb + iljnrgb
        if (imindist.lt.mindist) mindist = imindist
      enddo
      
c     apply cutoff if we have are doing molecule-molecule cutting
      if (cuttype.eq.MOLECULEMOLECULE) then
        if (mindist.lt.cut2) then
c         within cut distance (well, for this residue - we don't know
c         about other residues in the protein - this is a known BUG!)
          if (mindist.gt.feather2) then
c           however, the interaction must be scaled
            sclfac = featherfac * (cut2-mindist)
            nrg = nrg * sclfac
            cnrg = cnrg * sclfac
            ljnrg = ljnrg * sclfac
            cnrgf = cnrgf*sclfac
            cnrgb = cnrgb*sclfac
            ljnrgf = ljnrgf*sclfac
            ljnrgb = ljnrgb*sclfac
          endif
        else
c         cut the interaction
          nrg = ZERO
          cnrg = ZERO
          ljnrg = ZERO
          cnrgf = ZERO
          cnrgb = ZERO
          ljnrgf = ZERO
          ljnrgb = ZERO
        endif
      endif            

      residueProteinEnergy = nrg
      
      return
      
      end
