      double precision function soluteSingleSolventEnergy(isol,isvn,
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
c     between solute isol and solvent isvn
c
c     (C) Christopher Woods     29 January 2002
c
c##########################################################

      integer isol,isvn
      double precision cnrg,ljnrg,nrg,mindist
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb
      double precision distx,disty,distz,dist2,sclfac
      integer i,j,nat1,nat2,istart
      include 'nbenergyvar2.inc'
      double precision charge,sigma,epslon      
      double precision scl,sclf,sclb
      double precision getSolventRadius
        
      nrg = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      ljnrgf = ZERO
      cnrgb = ZERO
      ljnrgb = ZERO
      mindist = 999999.0
      soluteSingleSolventEnergy = ZERO
      softcore = .false.
      ddd = .false.
      if (IgnoreSolute(isol)) return

      nat1 = NSolAtoms(isol)
      nat2 = NSvnAtoms(isvn)
      
      if (nat1.le.0 .or. nat2.le.0) return
      
c     return if either the solvent or solute has been done before
c     (prevents double adding energies)
      if (DoneSolNrg(isol) .or. DoneSvnNrg(isvn)) return
      
c     now check the distance between the centre of the solvent, and centre
c     of the solute - also use this distance to select
c     the nearest image
      distx = dabs(SolCentre(isol,1)-SvnCentre(isvn,1))
      disty = dabs(SolCentre(isol,2)-SvnCentre(isvn,2)) 
      distz = dabs(SolCentre(isol,3)-SvnCentre(isvn,3))
            
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
      if (isol.eq.SvnPrefSol) then
        SvnNewProb(isvn) = dist2
      endif

c     cut the interaction if this distance is greater than the cutoff plus
c     the solute radius 
      if ( .not. SIRE_COMPATIBILITY_MODE ) then
        if (dist2 .gt. (cut+SolRadius(isol))**2) return
      else
        if (dist2 .gt. (cut+SolRadius(isol)+getSolventRadius(isvn))**2) return
      endif      
            
c     loop over all solvent-atom - solute-atom pairs, and calculate
c     thier interaction energy
      funcisol1 = isol
      funcisvn2 = isvn
c     Are we going to use a softcore ? 
      if (SoluteSoftCore(isol)) softcore = .true.
c     Activate Distance Dependent Dielectric calculation
      if (QuickMove .and. DDDon) ddd = .true. 
      do i=1,nat1
        funcat1 = i

        if ((SvnOptimization(isvn).eq.OPTTIP3P.or.SvnOptimization(isvn).eq.OPTTIP4P).and..not.doDerivative) then
c          if (isvn.eq.1) write(6,*)"Using optimized routines..."
c         Optimized for TIP4P-TIP4P interactions

c         LJ on oxygens only
          funcat2 = 1
          include 'nbinitial.inc'
          include 'nbsolsvn.inc'
          include 'nbenergy_lj.inc'         
        
          nrg = nrg + funcnbnrg
          ljnrg = ljnrg + funcljnrg
          ljnrgf = ljnrgf + funcljnrgf
          ljnrgb = ljnrgb + funcljnrgb
          if (funcmindist.lt.mindist) mindist = funcmindist

c         Charges on the rest of the particles, i.e., between hydrogens and extra point particle
          if (SvnOptimiZation(isvn).eq.OPTTIP3P) then
            istart = 1
          else
            istart = 2
          endif
          do j=istart,nat2
            funcat2 = j

            include 'nbinitial.inc'
            include 'nbsolsvn.inc'
            include 'nbenergy_coul.inc'
          
            nrg = nrg + funcnbnrg
            cnrg = cnrg + funccnrg
            cnrgf = cnrgf + funccnrgf
            cnrgb = cnrgb + funccnrgb

            if (funcmindist.lt.mindist) mindist = funcmindist
          enddo
        else
          do j=1,nat2
            funcat2 = j
  
            include 'nbinitial.inc'
            include 'nbsolsvn.inc'
            if (doDerivative) then
              include 'nbderivative.inc'
            else
              include 'nbenergy2.inc'
            endif

            nrg = nrg + funcnbnrg
            cnrg = cnrg + funccnrg
            ljnrg = ljnrg + funcljnrg
            cnrgf = cnrgf + funccnrgf
            cnrgb = cnrgb + funccnrgb
            ljnrgf = ljnrgf + funcljnrgf
            ljnrgb = ljnrgb + funcljnrgb
            if (funcmindist.lt.mindist) mindist = funcmindist
          enddo
        endif
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
      
      soluteSingleSolventEnergy = nrg
      
      return
      end
