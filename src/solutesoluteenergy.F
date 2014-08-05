      double precision function soluteSoluteEnergy(isol1,isol2,
     .           cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
      include 'dimensions.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'simparams.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'quickpot.inc'
c##############################################################
c
c     This function returns the total solute-solute interaction
c     energy between solutes isol1 and isol2
c
c     (C) Christopher Woods       29 January 2002
c
c##############################################################

      integer isol1,isol2
      double precision cnrg,ljnrg,nrg,mindist
      double precision cnrgf,ljnrgf,cnrgb,ljnrgb
      double precision distx,disty,distz,dist2,sclfac
      double precision scl,sclf,sclb
      double precision getLowScale,getHighScale
      integer i,j,nat1,nat2
      include 'nbenergyvar2.inc'
      double precision charge,sigma,epslon      

      soluteSoluteEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgb = ZERO
      cnrgf = ZERO
      ljnrgb = ZERO
      ljnrgf = ZERO
      nrg = ZERO
      mindist = 999999.0
      softcore = .false.
      ddd = .false.
      if (NSolutes.le.0  .or. isol1.eq.isol2 .or. isol1.gt.NSolutes
     .          .or. isol2.gt.NSolutes) then
        return
      else if (IgnoreSolute(isol1) .or. IgnoreSolute(isol2)) then
        return
c     return if either of the solutes has been done before (prevents double
c     calculation of shared energies)
      else if (DoneSolNrg(isol1) .or. DoneSolNrg(isol2)) then
        return
      else if (DualTopologySolute(isol1).eq.isol2 .or. DualTopologySolute(isol2).eq.isol1) then
c       do not calculate the energy between dual topology pairs
c        write(printstring,*) 'sol ',isol1,' and sol ',isol2,' NO NRG'
c        call printLine(DEBUG,printstring)
        return
      else if (DualTopologySolute(isol1).eq.0 .and. DualTopologySolute(isol2).eq.0) then
c       there are no dual topology solutes
        scl = ONE
        sclb = ONE
        sclf = ONE
c     there are dual topology solutes, so get the scale factors
      else if (DualTopologySolute(isol1).gt.0) then
c       solute 1 is a low dual topology solute
        scl = getLowScale(Lambda)
        sclb = getLowScale(LambdaB)
        sclf = getLowScale(LambdaF)
        doLowDerivative = .True.
        if (DualTopologySolute(isol2).gt.0) then
c         solute 2 is also a low dual topology solute
          scl = scl*scl
          sclb = sclb*sclb
          sclf = sclf*sclf
        else if (DualTopologySolute(isol2).lt.0) then
c         solute 2 is a high dual topology solute
          scl = scl * getHighScale(Lambda)
          sclb = sclb * getHighScale(LambdaB)
          sclf = sclf * getHighScale(LambdaF)
        endif
      else if (DualTopologySolute(isol1).eq.0) then
c       solute 1 is not a dual topology solute
        if (DualTopologySolute(isol2).gt.0) then
c         solute 2 is a low dual topology solute
          scl = getLowScale(Lambda)
          sclb = getLowScale(LambdaB)
          sclf = getLowScale(LambdaF)
          doLowDerivative = .True.
        else
c         solute 2 is a high dual topology solute
          scl = getHighScale(Lambda)
          sclb = getHighScale(LambdaB)
          sclf = getHighScale(LambdaF)
          doLowDerivative = .False.
        endif
      else
c       solute 1 is a high dual topology solute
        scl = getHighScale(Lambda)
        sclb = getHighScale(LambdaB)
        sclf = getHighScale(LambdaF)
        doLowDerivative = .False.
        if (DualTopologySolute(isol2).gt.0) then
c         solute 2 is a low dual topology solute
          scl = scl * getLowScale(Lambda)
          sclb = sclb * getLowScale(LambdaB)
          sclf = sclf * getLowScale(LambdaF)
        else if (DualTopologySolute(isol2).lt.0) then
c         solute 2 is also a high dual topology solute
          scl = scl*scl
          sclb = sclb*sclb
          sclf = sclf*sclf
        endif
      endif
          
      nat1 = NSolAtoms(isol1)
      nat2 = NSolAtoms(isol2)
      if (nat1.le.0 .or. nat2.le.0) return

c     check the distance between the solute centres
      distx = dabs(SolCentre(isol1,1)-SolCentre(isol2,1))
      disty = dabs(SolCentre(isol1,2)-SolCentre(isol2,2))
      distz = dabs(SolCentre(isol1,3)-SolCentre(isol2,3))
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
c     if this distance is greater than the cutoff, plus the solute radii
c     the cut the interaction      
      if (dist2 .gt. (cut+SolRadius(isol1)+SolRadius(isol2))**2) return

c     Are we going to use a softcore ? 
      if (SoluteSoftCore(isol1) .or. SoluteSoftCore(isol2)) softcore = .true.
      if (QuickMove .and. DDDon ) ddd = .true.
      funcisol1 = isol1
      funcisol2 = isol2
      do i=1,nat1
        funcat1 = i
        do j=1,nat2
          funcat2 = j
          
c         initialise the calculation
          include 'nbinitial.inc'
c         set up for a solute-solute interaction
          include 'nbsolsol.inc'
c         actually perform the calculation
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
      enddo

      if (cuttype.eq.RESIDUEMOLECULE .or. cuttype.eq.MOLECULEMOLECULE) then
        if (mindist.lt.cut2) then
c         within cutoff distance
          if (mindist.gt.feather2) then
c           in the feather region
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

c      write(printstring,*) 'isol1 ',isol1,' isol2 ',isol2,' scl ',scl,
c     .  ' cnrg ',cnrg,' cnrgf ',cnrgf,' cnrgb ',cnrgb
c      call printLine(DEBUG,printstring)
c      call printLine(WARNING,printstring)
c      if (QuickMove) then 
c         sclf = scl
c         sclb = scl
c      endif
c     scale the energies, unless we are calculating the lambda-gradient in case
c     the energies are already scaled
      if (.not.doDerivative) then
        nrg = nrg * scl
        cnrg = cnrg * scl
        ljnrg = ljnrg * scl
        cnrgf = cnrgf * sclf
        ljnrgf = ljnrgf * sclf
        cnrgb = cnrgb * sclb
        ljnrgb = ljnrgb * sclb
c       If softcore, the coulombic if scaled further
        if (Solutesoftcore(isol1) .or. SoluteSoftCore(isol2)) then
          cnrg = cnrg * (scl)**CoulSoftPower
          cnrgf = cnrgf * (sclf)**CoulSoftPower
          cnrgb = cnrgb * (sclb)**CoulSoftPower
        endif 
      endif      
c     soluteSoluteEnergy = nrg
      soluteSoluteEnergy = cnrg + ljnrg

      return
      end
