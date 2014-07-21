      double precision function GCsoluteSingleSolventEnergy(igcs,isvn,
     .                    cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist,
     .                    scl,sclf,sclb)
      include 'dimensions.inc'
      include 'simparams.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'quickpot.inc'
c##########################################################
c
c     This function calculates the interaction energy
c     between GCsolute igcs and solvent isvn
c
c     (C) MSB 22nd April 2009 
c
c##########################################################

      integer igcs,isvn
      double precision cnrg,ljnrg,nrg,mindist
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb
      double precision distx,disty,distz,dist2,sclfac
      integer i,j,nat1,nat2
      include 'nbenergyvar2.inc'
      double precision charge,sigma,epslon      
      double precision scl,sclf,sclb
      nrg = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      ljnrgf = ZERO
      cnrgb = ZERO
      ljnrgb = ZERO
      mindist = 999999.0
      GCsoluteSingleSolventEnergy = ZERO
      softcore = .false.
      ddd = .false.
c      if (IgnoreSolute(isol)) return

      nat1 = NGCSAtoms(igcs)
      nat2 = NSvnAtoms(isvn)
      
      if (nat1.le.0 .or. nat2.le.0) return
      
c     return if either the solvent or solute has been done before
c     (prevents double adding energies)
      if (DoneGCSNrg(igcs) .or. DoneSvnNrg(isvn)) return
      
c     now check the distance between the centre of the solvent, and centre
c     of the GC-solute - also use this distance to select
c     the nearest image
      distx = dabs(GCSCentre(igcs,1)-SvnCentre(isvn,1))
      disty = dabs(GCSCentre(igcs,2)-SvnCentre(isvn,2)) 
      distz = dabs(GCSCentre(igcs,3)-SvnCentre(isvn,3))
            
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

c     if the solvent molecule is interacting with the preferential
c     sampling solute, then calculate SvnNewProb for this solvent
c     molecule
c      if (isol.eq.SvnPrefSol) then
c        SvnNewProb(isvn) = dist2
c      endif
       if (dist2.gt.cut2) return
       

c     cut the interaction if this distance is greater than the cutoff plus
c     the solute radius 
c      if (dist2 .gt. (cut+GCSRadius(igcs))**2) return
                  
c     loop over all solvent-atom - solute-atom pairs, and calculate
c     thier interaction energy
      funcigcs1 = igcs
      funcisvn2 = isvn
c     Are we going to use a softcore ? 
      if (GCsoluteSoftCore(igcs)) softcore = .true.
c     Activate Distance Dependent Dielectric calculation
      if (QuickMove .and. DDDon) ddd = .true. 
      do i=1,nat1
        funcat1 = i

        do j=1,nat2
          funcat2 = j

          include 'nbinitial.inc'
          include 'nbgcssvn.inc'
          include 'nbenergy3.inc'

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

      if (cuttype.eq.RESIDUEMOLECULE .or. cuttype.eq.MOLECULEMOLECULE) then
        if (mindist.lt.cut2) then
c         smallest distance within cutoff
          if (mindist.gt.feather2) then
c           however, outside into feather region
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
          nrg = ZERO
          cnrg = ZERO
          ljnrg = ZERO
          cnrgf = ZERO
          ljnrgf = ZERO
          cnrgb = ZERO
          ljnrgb = ZERO
        endif
      endif

c       Now scale these energies according to theta...
            nrg = nrg * thetaparam(igcs)
            cnrg = cnrg * thetaparam(igcs)
            ljnrg = ljnrg * thetaparam(igcs)
            cnrgf = cnrgf * thetaparam(igcs)
            cnrgb = cnrgb * thetaparam(igcs)
            ljnrgf = ljnrgf * thetaparam(igcs)
            ljnrgb = ljnrgb * thetaparam(igcs)           
      
      GCsoluteSingleSolventEnergy = nrg
c      PRINT*, "mindist-nrg" ,mindist,GCsoluteSingleSolventEnergy
      return
      end
