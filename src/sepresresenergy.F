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

      sepResResEnergy = nrg
      
      return
      
      end
