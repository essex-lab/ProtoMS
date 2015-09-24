      double precision function getSolventRadius(isvn)
      include 'dimensions.inc'
      include 'solvents.inc'

      integer i,isvn
      double precision dx2, dy2, dz2
      double precision maxx2, maxy2, maxz2

      dx2 = 0
      dy2 = 0
      dz2 = 0
      maxx2 = 0
      maxy2 = 0
      maxz2 = 0

      do i=1,NSvnAtoms(isvn)
        dx2 = ( ZSvnCoords(isvn,i,1) - SvnCentre(isvn,1) )**2
        dy2 = ( ZSvnCoords(isvn,i,2) - SvnCentre(isvn,2) )**2
        dz2 = ( ZSvnCoords(isvn,i,3) - SvnCentre(isvn,3) )**2

        if (dx2 .gt. maxx2) maxx2 = dx2
        if (dy2 .gt. maxy2) maxy2 = dy2
        if (dz2 .gt. maxz2) maxz2 = dz2
      enddo

      getSolventRadius = sqrt(maxx2 + maxy2 + maxz2)

      return

      end

      double precision function solventSolventSingleEnergy(isvn1,isvn2,
     .           cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
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
c     This calculates the interaction energy between single solvent molecules
c     isvn1 and isvn2
c
c     (C) Christopher Woods   29 January 2002
c
c####################################################################

      integer isvn1,isvn2
      integer i,j,nats1,nats2,istart
      double precision mindist,dist2
      double precision distx,disty,distz
      double precision cnrg,ljnrg,nrg,sclfac
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb
      double precision getSolventRadius,svnradius1,svnradius2

      include 'nbenergyvar2.inc'
      double precision charge,sigma,epslon      
c     Set to one because there is no possible softcore
c     between two solvents at the current moment
      double precision scl,sclf,sclb
      scl = 1.0d0
      sclf = 1.0d0
      sclb = 1.0d0
c     end of cheat
      solventSolventSingleEnergy = ZERO
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
      if (DoneSvnNrg(isvn1) .or. DoneSvnNrg(isvn2)) return

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

      if (QuickMove.and. DDDon) ddd = .true.

      funcisvn1 = isvn1
      funcisvn2 = isvn2

      if ((SvnOptimization(isvn1).eq.OPTTIP4P.and.SvnOptimization(isvn1).eq.OPTTIP4P).or.
     .   (SvnOptimization(isvn1).eq.OPTTIP3P.and.SvnOptimization(isvn1).eq.OPTTIP3P)) then

c       Optimized for TIP4P-TIP4P interactions
c        write(6,*)"Using optimized routines"
c       LJ between oxygens only
        funcat1 = 1
        funcat2 = 1
        include 'nbinitial.inc'
        include 'nbsvnsvn.inc'
        include 'nbenergy_lj.inc'                 
        nrg = nrg + funcnbnrg
        ljnrg = ljnrg + funcljnrg
        ljnrgf = ljnrgf + funcljnrgf
        ljnrgb = ljnrgb + funcljnrgb
        if (SIRE_COMPATIBILITY_MODE .and. (funcmindist .lt. mindist)) 
     .          mindist = funcmindist

c       Coulomb between the rest of the particles, i.e., between hydrogens and extra point particle
        if (SvnOptimiZation(isvn1).eq.OPTTIP3P) then
          istart = 1
        else
          istart = 2
        endif
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(GUIDED)
!$OMP&PRIVATE(j,funcat1,funcat2,funcnbnrg,funccnrg,funcljnrg,funccnrgb,funcljnrgf,funcljnrgb)
!$OMP&PRIVATE(invfuncdist12b,invfuncdist12f,invfuncdist6b,invfuncdist6f,softdistb,softdistf)
!$OMP&PRIVATE(invfuncdistb,funcdistb6,invfuncdistf,funcdistf6,funcdistb2,funcdistf2)
!$OMP&PRIVATE(funcdzb,funcdzf,funcdyb,funcdyf,funcdxb,funcdxf,invfuncdist12,invfuncdist6)
!$OMP&PRIVATE(sigat2,sigat1,softpowerinv,invfuncdist,funcdist6,funcdist2)
!$OMP&PRIVATE(funcsig3b,funcsqrtepsb,funcsig3f,funcsqrtepsf)
!$OMP&PRIVATE(funcieps2f,funcieps2b,funcisig2f,funcisig2b,funcichg2b,funcichg2f)
!$OMP&PRIVATE(funcieps1f,funcieps1b,funcisig1f,funcisig1b,funcichg1f,funcichg1b)
!$OMP&PRIVATE(funcsig3,funcsqrteps)
!$OMP&PRIVATE(funcieps2,funcisig2,funcichg2,funcieps1,funcisig1)
!$OMP&PRIVATE(funcichg1,funcsigsclb,funcsigsclf,funcdogeom,funcisvn2,funcisvn1)
!$OMP&PRIVATE(funcisol2,funcisol1,funcdz,funcdy,funcdx,pertchg,shiftb,shiftf)
!$OMP&PRIVATE(shift,softdist,funcbijb,funcbijf,funcbij,funcaijf,funcaijb,funcaij)
!$OMP&PRIVATE(qiqjb,qiqjf,qiqj,ipar2,ipar1,funccnrgf,funcnbnrgb,funcnbnrgf,funcmindist)

!$OMP&SHARED(funciprot1,funciprot2,istart,nats2,ddd,softcore,scl,sclf,sclb,DeltaCore)
!$OMP&SHARED(DeltaCoreCoul,SoftType,boundary,boundx,boundy,boundz,boxdimension)
!$OMP&SHARED(AIJ,BIJ,CLJPair,ChgIJ,lambda,lambdaf,lambdab,DifferentLambdas)
!$OMP&SHARED(SIRE_COMPATIBILITY_MODE,lorentzberthelot,VariableEps,VariableSig,VariableChg)
!$OMP&SHARED(ZSvnPar,ZSvnCoords)

!$OMP&REDUCTION(+:nrg,cnrg,cnrgf,cnrgb)
!$OMP&REDUCTION(min: mindist)
        do i=istart,nats1
          funcat1 = i
          do j=istart,nats2
            funcat2 = j

            include 'nbinitial.inc'
            include 'nbsvnsvn.inc'
            include 'nbenergy_coul.inc'
          
            nrg = nrg + funcnbnrg
            cnrg = cnrg + funccnrg
            cnrgf = cnrgf + funccnrgf
            cnrgb = cnrgb + funccnrgb

            if (SIRE_COMPATIBILITY_MODE .and. (funcmindist .lt. mindist)) 
     .          mindist = funcmindist
          enddo
        enddo      
!$OMP END PARALLEL DO
 
      else

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(GUIDED)
!$OMP&PRIVATE(j,funcat1,funcat2,funcnbnrg,funccnrg,funcljnrg,funccnrgb,funcljnrgf,funcljnrgb)
!$OMP&PRIVATE(invfuncdist12b,invfuncdist12f,invfuncdist6b,invfuncdist6f,softdistb,softdistf)
!$OMP&PRIVATE(invfuncdistb,funcdistb6,invfuncdistf,funcdistf6,funcdistb2,funcdistf2)
!$OMP&PRIVATE(funcdzb,funcdzf,funcdyb,funcdyf,funcdxb,funcdxf,invfuncdist12,invfuncdist6)
!$OMP&PRIVATE(sigat2,sigat1,softpowerinv,invfuncdist,funcdist6,funcdist2)
!$OMP&PRIVATE(funcsig3b,funcsqrtepsb,funcsig3f,funcsqrtepsf)
!$OMP&PRIVATE(funcieps2f,funcieps2b,funcisig2f,funcisig2b,funcichg2b,funcichg2f)
!$OMP&PRIVATE(funcieps1f,funcieps1b,funcisig1f,funcisig1b,funcichg1f,funcichg1b)
!$OMP&PRIVATE(funcsig3,funcsqrteps)
!$OMP&PRIVATE(funcieps2,funcisig2,funcichg2,funcieps1,funcisig1)
!$OMP&PRIVATE(funcichg1,funcsigsclb,funcsigsclf,funcdogeom,funcisvn2,funcisvn1)
!$OMP&PRIVATE(funcisol2,funcisol1,funcdz,funcdy,funcdx,pertchg,shiftb,shiftf)
!$OMP&PRIVATE(shift,softdist,funcbijb,funcbijf,funcbij,funcaijf,funcaijb,funcaij)
!$OMP&PRIVATE(qiqjb,qiqjf,qiqj,ipar2,ipar1,funccnrgf,funcnbnrgb,funcnbnrgf,funcmindist)

!$OMP&SHARED(funciprot1,funciprot2,nats2,ddd,softcore,scl,sclf,sclb,DeltaCore)
!$OMP&SHARED(DeltaCoreCoul,SoftType,boundary,boundx,boundy,boundz,boxdimension)
!$OMP&SHARED(AIJ,BIJ,CLJPair,ChgIJ,lambda,lambdaf,lambdab,DifferentLambdas)
!$OMP&SHARED(SIRE_COMPATIBILITY_MODE,lorentzberthelot,VariableEps,VariableSig,VariableChg)
!$OMP&SHARED(ZSvnPar,ZSvnCoords)

!$OMP&REDUCTION(+:nrg,cnrg,ljnrg,cnrgf,cnrgb,ljnrgf,ljnrgb)
!$OMP&REDUCTION(min: mindist)
        do i=1,nats1
          funcat1 = i
          do j=1,nats2
            funcat2 = j

c           intialise the nb calculation          
            include 'nbinitial.inc'
c           set up for a solvent solvent interaction
            include 'nbsvnsvn.inc'
c           actually calculate the energy
            include 'nbenergy2.inc'
          
            nrg = nrg + funcnbnrg
            cnrg = cnrg + funccnrg
            ljnrg = ljnrg + funcljnrg
            cnrgf = cnrgf + funccnrgf
            cnrgb = cnrgb + funccnrgb
            ljnrgf = ljnrgf + funcljnrgf
            ljnrgb = ljnrgb + funcljnrgb
            if (SIRE_COMPATIBILITY_MODE .and. (funcmindist .lt. mindist)) 
     .          mindist = funcmindist
          enddo
        enddo
!$OMP END PARALLEL DO
      endif

      if (cuttype.eq.RESIDUEMOLECULE .or. cuttype.eq.MOLECULEMOLECULE) then
        if (mindist.gt.cut2) then
          nrg = 0
          cnrg = 0
          ljnrg = 0
          cnrgf = 0
          cnrgb = 0
          ljnrgf = 0
          ljnrgb = 0
        
        else if (mindist.gt.feather2) then
          
c          print *,sqrt(mindist)," ",cnrg," ",ljnrg
c          print *,featherfac * (cut2-mindist)

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
     
      solventSolventSingleEnergy = nrg
      
      return
      end
