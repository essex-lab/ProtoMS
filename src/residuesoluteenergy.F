      double precision function residueSoluteEnergy(iprot,ires,isol,
     .           cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist,
     .           scl,sclf,sclb)
      include 'dimensions.inc'
      include 'proteins.inc'
      include 'simparams.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'quickpot.inc'
c#######################################################################
c
c     This returns the interaction energy between residue ires of 
c     protein iprot, with solute isol
c
c     (C) Christopher Woods     29 January 2002
c
c#######################################################################

      integer iprot,ires,isol
      double precision cnrg,ljnrg,nrg
      double precision cnrgf,ljnrgf,cnrgb,ljnrgb
      double precision mindist,sclfac
      double precision distx,disty,distz,dist2
      integer i,j,nat,resstrt,resend
      include 'nbenergyvar2.inc'
      double precision charge,sigma,epslon      
      double precision scl,sclf,sclb
      residueSoluteEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      nrg = ZERO
      mindist = 90000.0
      softcore = .false.
      ddd = .false.
      if (NProteins.le.0 .or. iprot.gt.NProteins .or. NSolutes.le.0
     .       .or. isol.gt.NSolutes) return
                   
      if (IgnoreSolute(isol)) return
     
      nat = NSolAtoms(isol)
      if (nat.le.0) return

c     return if either the solute or residue has been done before (this
c     prevents double adding energies)
      if (DoneResNrg(iprot,ires) .or. DoneSolNrg(isol)) return

c     check the distance between the centre of the residue and the centre
c     of the solute
      distx = dabs(ResCentre(iprot,ires,1)-SolCentre(isol,1))
      disty = dabs(ResCentre(iprot,ires,2)-SolCentre(isol,2))
      distz = dabs(ResCentre(iprot,ires,3)-SolCentre(isol,3))
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

c     if the distance is greater than cut+solRadius+resRadius then cut
      if (dist2 .gt. (cut+SolRadius(isol)+ResRadius(iprot,ires)+ONE)**2) return

c     assume residue runs concurrently from N to O
      resstrt = BBZid(iprot,ires,1)
      resend = BBZid(iprot,ires,4)
      
      funciprot1 = iprot
      funcisol2 = isol
c     Are we going to use a softcore ? 
      if (SoluteSoftCore(isol)) softcore = .true.
c     Activate Distance Dependent Dielectric calculation
      if (QuickMove .and. DDDon) ddd = .true.
      do i=resstrt,resend
        funcat1 = i  
        do j=1,nat
          funcat2 = j        

          include 'nbinitial.inc'
          include 'nbprosol.inc'
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

      if (cuttype.eq.RESIDUEMOLECULE) then
c     performing protein-residue - solute molecule cutting...
        if (mindist.lt.cut2) then
          if (mindist.gt.feather2) then
            sclfac = featherfac * (cut2-mindist)
            nrg = nrg * sclfac
            cnrg = cnrg * sclfac
            ljnrg = ljnrg * sclfac
            cnrgf = cnrgf * sclfac
            ljnrgf = ljnrgf * sclfac
            cnrgb = cnrgb * sclfac
            ljnrgb = ljnrgb * sclfac
          endif
        else
          nrg = ZERO
          cnrg = ZERO
          ljnrg = ZERO
          cnrgf = ZERO
          cnrgb = ZERO
          ljnrgf = ZERO
          ljnrgb = ZERO
        endif
      endif
      
      residueSoluteEnergy = nrg

      return
      end
