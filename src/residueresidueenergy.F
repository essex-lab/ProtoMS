      double precision function residueResidueEnergy(iprot,ires1,ires2,
     .                   cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'flags.inc'
c####################################################################
c
c     This function calculates the intramolecular protein
c     residue-residue non-bonded energy...
c
c     It also returns the coulomb and lj components in cnrg and ljnrg
c
c     (C) Christopher Woods  25 January 2002
c
c####################################################################

      integer iprot,ires1,ires2
      double precision cnrg,ljnrg,nrg,cnrgf,cnrgb,ljnrgf,ljnrgb,mindist
      double precision sepResResEnergy,bondedResResEnergy
      double precision sclfac,dist2,distx,disty,distz            
      include 'nbenergyvar2.inc'
      
      residueResidueEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      ljnrgf = ZERO
      cnrgb = ZERO
      ljnrgb = ZERO
      mindist = 999990.0
      softcore = .false.
      if (ires1.eq.ires2) return

c     if we are doing residue-molecule cutting, then check the distance
c     between the centres of the two residues
      if (cuttype.eq.RESIDUEMOLECULE) then
        distx = dabs(ResCentre(iprot,ires1,1)-ResCentre(iprot,ires2,1))
        disty = dabs(ResCentre(iprot,ires1,2)-ResCentre(iprot,ires2,2))
        distz = dabs(ResCentre(iprot,ires1,3)-ResCentre(iprot,ires2,3))
        if (boundary.eq.PERIODICBOUNDARIES) then
c         since this is the *intramolecular* energy, it should not be subjected to
c         periodic boundaries - thus boundx/boundy/boundz should be false
          boundx = .false.
          boundy = .false.
          boundz = .false.
        endif
        dist2 = distx**2 + disty**2 + distz**2
c       now if the distance plus the radii of the two residues is greater than
c       the cutoff, make a cut!
        mindist = dist2
        if (dist2 .gt. (cut+ResRadius(iprot,ires1)+ResRadius(iprot,ires2)+ONE)**2) return
      endif
      
c     if the energy of one of the residues has already been done, then
c     return (this prevents double adding energies!
      if (DoneResNrg(iprot,ires1) .or. DoneResNrg(iprot,ires2)) return
      
c     check to see if the residues are bonded...
      if (ResResBond(iprot,ires1,ires2)) then
c       they are bonded...
        nrg = bondedResResEnergy(iprot,ires1,ires2,cnrg,ljnrg,cnrgb,ljnrgb,
     .                                                    cnrgf,ljnrgf,mindist)
      else 
c       they are not bonded
        nrg = sepResResEnergy(iprot,ires1,ires2,cnrg,ljnrg,cnrgb,ljnrgb,
     .                                                     cnrgf,ljnrgf,mindist)
      endif
      
      if (cuttype.eq.RESIDUEMOLECULE) then
c       doing protein-residue - protein-residue cutting
        if (mindist.lt.cut2) then
c         smallest distance within cutoff...
          if (mindist.gt.feather2) then
c           however, it is greater that the feather distance...
c           thus the interaction must be scaled (feathered)
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
c         falls outside cutoff, so must be set to zero...
          nrg = ZERO
          cnrg = ZERO
          ljnrg = ZERO
          cnrgb = ZERO
          cnrgf = ZERO
          ljnrgf = ZERO
          ljnrgb = ZERO
        endif
      endif

      residueResidueEnergy = nrg
      return
      
      end
