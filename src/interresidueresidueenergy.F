      double precision function interResidueResidueEnergy(iprot1,ires1,iprot2,ires2,
     .                             cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'proteins.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'quickpot.inc'
c#######################################################
c
c     This function calculates the inter-protein
c     residue-residue energy
c
c     (C) Christopher Woods   25 January 2002
c
c########################################################

      integer iprot1,iprot2,ires1,ires2,i,j
      integer resstrt1,resstrt2,resend1,resend2
      double precision nrg,cnrg,ljnrg,sclfac,cnrgf,cnrgb,ljnrgf,ljnrgb,mindist
      double precision residueResidueEnergy,dist2,distx,disty,distz
      include 'nbenergyvar2.inc'
      double precision charge,sigma,epslon      
c     temporary cheat
      double precision scl,sclf,sclb
      scl = 1.0d0
      sclf = 1.0d0
      sclb = 1.0d0
c     end of cheat
      interResidueResidueEnergy = ZERO
      nrg = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      mindist = 99999.0
      softcore = .false.
      ddd = .false.
      if (iprot1.eq.iprot2) then
        interResidueResidueEnergy = residueResidueEnergy(iprot1,ires1,ires2,cnrg,ljnrg,
     .                                    cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
        return
      endif
      
      if (iprot1.le.0 .or. iprot1.gt.NProteins .or. iprot2.le.0 .or. iprot2.gt.NProteins) return
      if (ires1.le.0 .or. ires1.gt.NResidues(iprot1)
     .     .or. ires2.le.0 .or. ires2.gt.NResidues(iprot2)) return
     
      if (DoneResNrg(iprot1,ires1) .or. DoneResNrg(iprot2,ires2)) return
      
c     use the distance between the center of the residues to decide whether
c     to use periodic boundaries      
      distx = dabs(ResCentre(iprot1,ires1,1)-ResCentre(iprot2,ires2,1))
      disty = dabs(ResCentre(iprot1,ires1,2)-ResCentre(iprot2,ires2,2))
      distz = dabs(ResCentre(iprot1,ires1,3)-ResCentre(iprot2,ires2,3))
      if (boundary.eq.PERIODICBOUNDARIES) then
        boundx = .false.
        boundy = .false.
        boundz = .false.
        if (distx.gt.BoxDimension2(1)) then
          distx = distx - BoxDimension(1)
          boundx = .true.
        endif
        if (disty.gt.BoxDimension2(2)) then
          disty = disty - BoxDimension(2)
          boundy = .true.
        endif
        if (distz.gt.BoxDimension2(3)) then
          distz = distz - BoxDimension(3)
          boundz = .true.
        endif
      endif
      dist2 = distx**2 + disty**2 + distz**2
      mindist = dist2
      if (cuttype.eq.RESIDUEMOLECULE) then
c       now if the distance plus the radii of the two residues is greater than
c       the cutoff, make a cut!
        if (dist2 .gt. (cut+ResRadius(iprot1,ires1)+ResRadius(iprot2,ires2))**2) return
      endif        
     
c     assume residues run cuncurrently from N to O     
      resstrt1 = BBZid(iprot1,ires1,1)
      resend1 = BBZid(iprot1,ires1,4)
      resstrt2 = BBZid(iprot2,ires2,1)
      resend2 = BBZid(iprot2,ires2,4)

      if (QuickMove .and. DDDon) ddd = .true.
      
      funciprot1 = iprot1
      funciprot2 = iprot2
      
      do i=resstrt1,resend1
        funcat1 = i
        
        do j=resstrt2,resend2
          funcat2 = j

          include 'nbinitial.inc'
          include 'nbpropro.inc'
          include 'nbenergy2.inc'          

          nrg = nrg + funcnbnrg
          cnrg = cnrg + funccnrg
          ljnrg = ljnrg + funcljnrg
          cnrgf = cnrgf + funccnrgf
          cnrgb = cnrgb + funccnrgb
          ljnrgf = ljnrgf + funcljnrgf
          ljnrgb = ljnrgb + funcljnrgb
          if (funcmindist.lt.mindist) mindist = funcmindist
        enddo
      enddo

c     apply cutoffs if necessary...
      if (cuttype.eq.1) then
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
            cnrgf = cnrgf * sclfac
            cnrgb = cnrgb * sclfac
            ljnrgf = ljnrgf * sclfac
            ljnrgb = ljnrgb * sclfac
          endif
        else
c         falls outside cutoff, so must be set to zero...
          nrg = ZERO
          cnrg = ZERO
          ljnrg = ZERO
          cnrgf = ZERO
          cnrgb = ZERO
          ljnrgf = ZERO
          ljnrgb = ZERO
        endif
      endif
      
      interResidueResidueEnergy = nrg
      return
      end
