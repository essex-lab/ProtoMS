
      subroutine solventSolventEleEnergy(isvn1,isvn2,alpha,
     .           coulomb,realspace)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'enums.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'quickpot.inc'
c####################################################################
c
c     This calculates the electrostatic interaction energy between
c      single solvent molecules isvn1 and isvn2
c      both bare Coulomb and real-space part of PME is computed
c
c     (C) Samuel Genheden 27 September 2013
c
c####################################################################

      integer isvn1,isvn2
      integer i,j,nats1,nats2
      double precision mindist,dist2,alpha
      double precision distx,disty,distz,dist
      double precision coulomb,realspace
      double precision getSolventRadius,svnradius1,svnradius2,dx,dy,dz

      logical boundx,boundy,boundz
      double precision charge    
      
      coulomb = 0.0d0
      realspace = 0.0d0      

c     now check the distance between the *centres* of the solvent molecules
c     is less than cutoff      

c     also, solvent are handled slightly differently - unlike everything
c     else, which does periodic boundaries on an atom-atom basis, the solvent
c     does periodic boundaries based on a single point, in this case the centre
      distx = dabs(SvnCentre(isvn1,1)-SvnCentre(isvn2,1)) 
      disty = dabs(SvnCentre(isvn1,2)-SvnCentre(isvn2,2))
      distz = dabs(SvnCentre(isvn1,3)-SvnCentre(isvn2,3)) 
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

      if ( .not. SIRE_COMPATIBILITY_MODE ) then
        if ( mindist .gt. cut2 ) return
      else 
        svnradius1 = getSolventRadius(isvn1)
        svnradius2 = getSolventRadius(isvn2)

        if ( mindist .gt. (cut+svnradius1+svnradius2)**2 ) return
      endif 

      nats1 = NSvnAtoms(isvn1)
      nats2 = NSvnAtoms(isvn2)
      if (nats1.le.0  .or. nats2.le.0) return

      do i=1,nats1
        do j=1,nats2
          dx = dabs(ZSvnCoords(isvn1,i,1)-ZSvnCoords(isvn2,j,1))
          dy = dabs(ZSvnCoords(isvn1,i,2)-ZSvnCoords(isvn2,j,2))
          dz = dabs(ZSvnCoords(isvn1,i,3)-ZSvnCoords(isvn2,j,3))
          if (boundary.eq.PERIODICBOUNDARIES) then
            if (boundx) dx = dx - BoxDimension(1)
            if (boundy) dy = dy - BoxDimension(2)
            if (boundz) dz = dz - BoxDimension(3)         
          endif

          dist2 = dx**2+dy**2+dz**2
          dist = sqrt(dist2)
          realspace = realspace + charge(ZSvnPar(isvn1,i))*charge(ZSvnPar(isvn2,j))*erfc(alpha*dist)/dist         
          coulomb = coulomb + charge(ZSvnPar(isvn1,i))*charge(ZSvnPar(isvn2,j))/dist
        enddo
      enddo
     
      realspace = realspace*FAC4PIEPS0
      coulomb = coulomb*FAC4PIEPS0

      return
      end
