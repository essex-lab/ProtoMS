      double precision function residueSingleGCsoluteEnergy(iprot,ires,igcs,
     .             cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
      include 'dimensions.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'simparams.inc'
      include 'solvents.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'quickpot.inc'
      include 'solutes.inc'
c####################################################################
c
c     This calculates the interaction energy between residue ires
c     of protein iprot, and GCsolute molecule igcs
c
c     (C) MSB 22nd April 2009
c
c####################################################################

      integer iprot,ires,igcs
      double precision cnrg,ljnrg,nrg
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb
      double precision distx,disty,distz,dist2
      integer i,j,nat,resstrt,resend
      double precision sclfac
      double precision mindist
      include 'nbenergyvar2.inc'
      double precision charge,sigma,epslon      
c     temporary cheat

      residueSingleGCsoluteEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      nrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      mindist = 99900.0
      nat = NGCSAtoms(igcs)
      softcore = .false.
      ddd = .false.
      if (nat.le.0) return
        
c     return if either the residue or GC-solute has been done before
c     Might need to do something with this...?
c     COMMENTED OUT FOR THE MOMENT
      if (DoneResNrg(iprot,ires) .or. DoneGCSNrg(igcs)) return

c     now check the distance between the residue centre and centre
c     of the GC-solute molecule - cut if it greater than cutoff+residueRadius
      distx = dabs(ResCentre(iprot,ires,1)-GCSCentre(igcs,1))
      disty = dabs(ResCentre(iprot,ires,2)-GCSCentre(igcs,2))
      distz = dabs(ResCentre(iprot,ires,3)-GCSCentre(igcs,3))
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
c      PRINT*, "mindist" ,mindist
      if (dist2 .gt. (cut+ResRadius(iprot,ires))**2) return
c      PRINT*, "continue"
c     assume residue runs concurrently from N to O
      resstrt = BBZid(iprot,ires,1)
      resend = BBZid(iprot,ires,4)

      funciprot1 = iprot
      funcigcs2 = igcs
      if (QuickMove .and. DDDon) ddd = .true.
      if (GCsoluteSoftCore(igcs)) softcore = .true.
      do i=resstrt,resend
        funcat1 = i
              
        do j=1,nat
          funcat2 = j

          include 'nbinitial.inc'
          include 'nbprogcs.inc'
          include 'nbenergy3.inc'
          
          nrg = nrg + funcnbnrg
          cnrg = cnrg + funccnrg
c          PRINT*, "cnrg" ,cnrg
          ljnrg = ljnrg + funcljnrg
c           PRINT*, "ljnrg" ,ljnrg
          cnrgf = cnrgf + funccnrgf
          cnrgb = cnrgb + funccnrgb
          ljnrgf = ljnrgf + funcljnrgf
          ljnrgb = ljnrgb + funcljnrgb
          if (funcmindist.lt.mindist) mindist = funcmindist
        enddo
      enddo

      if (cuttype.eq.RESIDUEMOLECULE) then
c       doing protein-residue - GC-solute molecule cutting
        if (mindist.lt.cut2) then
          if (mindist.gt.feather2) then
c           must scale the interactions..
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
          sclfac = ZERO
        endif
      endif
      
c      Now scale the energies according to the value of theta...          
            nrg = nrg * thetaparam(igcs)
            cnrg = cnrg * thetaparam(igcs)
            ljnrg = ljnrg * thetaparam(igcs)            
            cnrgf = cnrgf * thetaparam(igcs)
            ljnrgf = ljnrgf * thetaparam(igcs)
            cnrgb = cnrgb * thetaparam(igcs)
            ljnrgb = ljnrgb * thetaparam(igcs)      
c      PRINT*, "cnrg" ,cnrg
      residueSingleGCsoluteEnergy = nrg
c      PRINT*,"residue-gc nrg" ,residueSingleGCsoluteEnergy
      return
      end
