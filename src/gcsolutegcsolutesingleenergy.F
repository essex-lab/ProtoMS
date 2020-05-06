      double precision function GCsoluteGCsoluteSingleEnergy(igcs1,igcs2,
     .           cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'solutes.inc'
      include 'enums.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'quickpot.inc'
      include 'stacks.inc'
c####################################################################
c
c     This calculates the interaction energy between single solute molecules
c     igcs1 and igcs2. In order to use this for the JAWS algorithm, there's a
c     dual theta scaling effect. Not sure if this will work, but we'll see!
c
c     (C) MSB 22nd April 09, modified 23rd October 2009
c
c####################################################################

      integer igcs1,igcs2
      integer i,j,nats1,nats2
      double precision mindist,dist2
      double precision distx,disty,distz
      double precision cnrg,ljnrg,nrg,sclfac
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb
      include 'nbenergyvar2.inc'
      double precision charge,sigma,epslon      

      GCsoluteGCsoluteSingleEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      nrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      mindist = 9999990.0
      softcore = .false.
      ddd = .false.
c     Might need a cheat or something here!      
c      if (DoneSvnNrg(igcs1)) print*, "START",igcs1

      if (DoneGCSNrg(igcs1) .and. DoneGCSNrg(igcs2)) return
c      if (DoneSvnNrg(igcs1)) print*, "START",igcs1
c      if (DoneSvnNrg(igcs2)) print*, "START",igcs2
    
c     now check the distance between the *centres* of the solute molecules
c     is less than cutoff      

c     also, solvent are handled slightly differently - unlike everything
c     else, which does periodic boundaries on an atom-atom basis, the solvent
c     does periodic boundaries based on a single point, in this case the centre
      distx = dabs(GCSCentre(igcs1,1)-GCSCentre(igcs2,1)) 
      disty = dabs(GCSCentre(igcs1,2)-GCSCentre(igcs2,2))
      distz = dabs(GCSCentre(igcs1,3)-GCSCentre(igcs2,3)) 
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
c      PRINT*, "dist2-cut2" ,dist2,cut2
      mindist = dist2
      if (dist2.gt.cut2) return
      
      nats1 = NGCSAtoms(igcs1)
      nats2 = NGCSAtoms(igcs2)
      if (nats1.le.0  .or. nats2.le.0) return

      if (QuickMove.and. DDDon) ddd = .true.
      if (GCsoluteSoftCore(igcs1) .or. GCsoluteSoftCore(igcs2)) softcore = .true.
      funcigcs1 = igcs1
      funcigcs2 = igcs2
      do i=1,nats1
        funcat1 = i
        do j=1,nats2
          funcat2 = j

c         intialise the nb calculation          
          include 'nbinitial.inc'
c         set up for a solute-solute interaction
          include 'nbgcsgcs.inc'
c         actually calculate the energy
          include 'nbenergy3.inc'
          
          nrg = nrg + funcnbnrg
          cnrg = cnrg + funccnrg
          ljnrg = ljnrg + funcljnrg
          cnrgf = cnrgf + funccnrgf
          cnrgb = cnrgb + funccnrgb
          ljnrgf = ljnrgf + funcljnrgf
          ljnrgb = ljnrgb + funcljnrgb
        enddo
      enddo

      if (cuttype.eq.RESIDUEMOLECULE .or. cuttype.eq.MOLECULEMOLECULE) then
c       also, unlike everything else, the cutting decision for the solvent
c       is made based on distance between the centre atoms - since we
c       have already eliminated the molecule if it cuts, we just need to
c       check if the interaction needs feathering
        if (mindist.gt.feather2) then
          sclfac = featherfac * (cut2-mindist)
          nrg = nrg * sclfac
          cnrg = cnrg * sclfac
          ljnrg = ljnrg * sclfac
          cnrgf = cnrgf*sclfac
          cnrgb = cnrgb*sclfac
          ljnrgf = ljnrgf*sclfac
          ljnrgb = ljnrgb*sclfac
        endif
      endif
c         print*, "pre-scaling" ,nrg
c      Finally, these energies need to be scaled according to theta1 and theta2
c      This is since the interaction energy needs to include BOTH of them...
          nrg = nrg * thetaparam(igcs1)*thetaparam(igcs2)
          cnrg = cnrg * thetaparam(igcs1)*thetaparam(igcs2)
          ljnrg = ljnrg * thetaparam(igcs1)*thetaparam(igcs2)
c          print*, "LJ", ljnrg
          cnrgf = cnrgf*thetaparam(igcs1)*thetaparam(igcs2)
          cnrgb = cnrgb*thetaparam(igcs1)*thetaparam(igcs2)
          ljnrgf = ljnrgf*thetaparam(igcs1)*thetaparam(igcs2)
          ljnrgb = ljnrgb*thetaparam(igcs1)*thetaparam(igcs2)
      

      if (doREST) then
        cnrg = cnrg * restscaling(5)
        ljnrg = ljnrg * restscaling(5)
        nrg = cnrg+ljnrg
        cnrgf = cnrgf * restscaling(5)
        ljnrgf = ljnrgf * restscaling(5)
        cnrgb = cnrgb * restscaling(5)
        ljnrgb = ljnrgb * restscaling(5)
      endif

c      PRINT*, "THETA-GCS", thetaparam(igcs1),thetaparam(igcs2)
      GCsoluteGCsoluteSingleEnergy = nrg
c      PRINT*, "Post-Scaling" ,igcs1,igcs2,GCsoluteGCsoluteSingleEnergy
      return
      end
