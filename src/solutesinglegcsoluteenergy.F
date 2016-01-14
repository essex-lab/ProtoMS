      double precision function soluteSingleGCSoluteEnergy(isol,igcs,
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
c ### delete later below
      include 'parallel.inc'
c##########################################################
c
c     This function calculates the interaction energy
c     between solute isol and gcsolute igcs
c
c     (C) MSB     8/12/09
c
c##########################################################

      integer isol,igcs
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
      soluteSingleGCSoluteEnergy = ZERO
      softcore = .false.
      ddd = .false.
      if (IgnoreSolute(isol)) return

      nat1 = NSolAtoms(isol)
      nat2 = NGCSAtoms(igcs)
      
      if (nat1.le.0 .or. nat2.le.0) return
      
c     return if either the solvent or solute has been done before
c     (prevents double adding energies)
      if (DoneSolNrg(isol) .or. DoneGCSNrg(igcs)) return
      
c     now check the distance between the centre of the solvent, and centre
c     of the solute - also use this distance to select
c     the nearest image
      distx = dabs(SolCentre(isol,1)-GCSCentre(igcs,1))
      disty = dabs(SolCentre(isol,2)-GCSCentre(igcs,2)) 
      distz = dabs(SolCentre(isol,3)-GCSCentre(igcs,3))
            
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
c     sampling solute, then calculate GCSNewProb for this solvent
c     molecule
c      if (isol.eq.GCSPrefSol) then
c        GCSNewProb(iGCS) = dist2
c      endif

c     cut the interaction if this distance is greater than the cutoff plus
c     the solute radius 
      if (dist2 .gt. (cut+SolRadius(isol))**2) return
                  
c     loop over all solvent-atom - solute-atom pairs, and calculate
c     thier interaction energy
      funcisol1 = isol
      funcigcs2 = igcs
c     Are we going to use a softcore ? 
      if (SoluteSoftCore(isol).and.GCsoluteSoftCore(igcs)) then
c       Break if both solute and gcsolute have softcore assigned
        call printLine(FATAL,"The calculation of interactions between softcored solute and softcored gcsolute is not implemented.")
        call printLine(FATAL,"The combination of lambda and JAWS calculations with softcores may be causing the issue.")
        call closeProgram(100,"")
      elseif  (SoluteSoftCore(isol)) then 
        softcore = .true.
        thet = scl
        thetf = sclf
        thetb = sclb
c        write(6,*)"scl, sclf, sclb: ",scl,sclf,sclb
c        thet = ONE-lambda
      elseif (GCsoluteSoftCore(igcs)) then
        softcore = .true.
        thet = thetaparam(funcigcs2)
        thetf = thet
        thetb = thet
      endif
c     Activate Distance Dependent Dielectric calculation
      if (QuickMove .and. DDDon) ddd = .true. 
      do i=1,nat1
        funcat1 = i

c      write(6,*) myid," Do we have softcore in solute-gcsolute interactions? ",softcore

        do j=1,nat2
          funcat2 = j

          include 'nbinitial.inc'
          include 'nbsolgcs.inc'
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
 
c        scale to the theta jaws parameter
            nrg = nrg * thetaparam(igcs)
            cnrg = cnrg * thetaparam(igcs)
            ljnrg = ljnrg * thetaparam(igcs)
            cnrgf = cnrgf * thetaparam(igcs)
            cnrgb = cnrgb * thetaparam(igcs)
            ljnrgf = ljnrgf * thetaparam(igcs)
            ljnrgb = ljnrgb * thetaparam(igcs)


      if (doREST) then
        cnrg = cnrg * restscaling(5)
        ljnrg = ljnrg * restscaling(5)
        nrg = cnrg + ljnrg
        cnrgf = cnrgf * restscaling(5)
        ljnrgf = ljnrgf * restscaling(5)
        cnrgb = cnrgb * restscaling(5)
        ljnrgb = ljnrgb * restscaling(5)
      endif
     
      soluteSingleGCSoluteEnergy = nrg

      if (isol.eq.1) then
        if (myid.eq.1) write(6,*) myid," gcsol-solute energy is ", soluteSingleGCSoluteEnergy,". cnrg ",cnrg,
     .                            " Mylambda ",lambdaval(mylambda)
      endif
      return
      end
