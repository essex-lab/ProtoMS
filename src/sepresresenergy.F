      double precision function sepResResEnergy(iprot,ires1,ires2,cnrg,ljnrg,
     .                                           cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'printer.inc'
      include 'enums.inc'
      include 'flags.inc'
      include 'parameters.inc'
      include 'simparams.inc'
      include 'quickpot.inc'
c#######################################################
c
c     This is an internal function only intended to be
c     used by residueResidueEnergy, used to calculate
c     the residue-residue energy between residues within
c     a protein that are not bonded...
c
c     (C) Christopher Woods   25 January 2002
c
c########################################################

      integer iprot,ires1,ires2
      double precision nrg,cnrg,ljnrg,cnrgb,cnrgf,ljnrgb,ljnrgf,mindist
      double precision getLowScale,scl,sclf,sclb
      integer zstrt1,zend1,zstrt2,zend2
      integer i,j
      include 'nbenergyvar2.inc'
      double precision charge,sigma,epslon      

      sepResResEnergy = ZERO
      nrg = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      mindist = 99999.0
      softcore = .false.
      scl = 1.0d0
      sclf = 1.0d0
      sclb = 1.0d0

      ddd = .false.

c     part of the promc file spec is that an amino-acid starts with
c     a N, and ends with the O. Since the zmatrix ID is sequential,
c     we can use that to run over the two residues...
      zstrt1 = BBZid(iprot,ires1,1)
      zend1 = BBZid(iprot,ires1,4)
      zstrt2 = BBZid(iprot,ires2,1)
      zend2 = BBZid(iprot,ires2,4)
      
      funciprot1 = iprot
      funciprot2 = iprot
      if (QuickMove .and. DDDon) ddd = .true.

c     are we softening the residue-residue interactions?
      if (ProteinSoftCore) then
        softcore = .true.
        scl = getLowScale(Lambda)
        sclf = getLowScale(LambdaF)
        sclb = getLowScale(LambdaB)
        print *,softcore,DeltaCore,scl,sclf,sclb

      endif

c     now run over all pairs...
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

!$OMP&SHARED(funciprot1,funciprot2,zstrt2,zend2,ddd,softcore,scl,sclf,sclb,DeltaCore)
!$OMP&SHARED(DeltaCoreCoul,SoftType,boundary,boundx,boundy,boundz,boxdimension)
!$OMP&SHARED(AIJ,BIJ,CLJPair,ChgIJ,lambda,lambdaf,lambdab,DifferentLambdas)
!$OMP&SHARED(SIRE_COMPATIBILITY_MODE,lorentzberthelot,VariableEps,VariableSig,VariableChg)
!$OMP&SHARED(ZCoords,ZPar)

!$OMP&REDUCTION(+:nrg, cnrg, ljnrg,cnrgf,cnrgb,ljnrgf,ljnrgb)
!$OMP&REDUCTION(min: mindist)
      do i=zstrt1,zend1
        funcat1 = i
        do j=zstrt2,zend2
          funcat2 = j

          include 'nbinitial.inc'
          include 'nbpropro.inc'
          include 'nbenergy2.inc'

c         output to funcnbnrg,funccnrg,funcljnrg,funcmindist
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
!$OMP END PARALLEL DO

10    format('Energies= ',8f8.3)
c     write(printstring,10) nrg,cnrg,ljnrg,cnrgf,cnrgb,ljnrgf,ljnrgb,mindist
c     call printLine(RESULTS,printstring)

      sepResResEnergy = nrg
      
      return
      
      end
