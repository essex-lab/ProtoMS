      double precision function sepResResEnergy(iprot,ires1,ires2,cnrg,ljnrg,
     .     cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
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
c     Extensive rewrite to improve optimisation of code
c
c     (C) Chris Cave-Ayland March 2015     
c
c########################################################

      integer iprot,ires1,ires2,iprot1,iprot2
      double precision nrg,cnrg,ljnrg,cnrgb,cnrgf,ljnrgb,ljnrgf,mindist
      double precision getLowScale,scl,sclf,sclb
      integer zstrt1,zend1,zstrt2,zend2
      integer i,j,ui,uj,tmpi,tmpj, ind, k, par
      include 'nbenergyvar2.inc'
      double precision charge,sigma,epslon
      double precision invr2, tmp2, tmp3, tmp
      double precision invrs(600),invrs6(600)

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

      if (ddd .or. ProteinSoftCore) then
c     run unoptimised version of nbenergy as this can
c     handle softcores and implicit solvent
c     now run over all pairs...
         do i=zstrt1,zend1
            funcat1 = i
            do j=zstrt2,zend2
               funcat2 = j
               
               include 'nbinitial.inc'
               include 'nbpropro.inc'
               include 'nbenergy2.inc'
               
c     output to funcnbnrg,funccnrg,funcljnrg,funcmindist
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
         
      else
c     optimised version of nbenergy cannot handle protein-protein
c     softcores or implicitt solvent
         ui = zend1-zstrt1+1
         uj = zend2-zstrt2+1

         tmpi = zstrt1
         tmpj = zstrt2

         ind=1
         do i=1,ui
            do j=1,uj
               invrs(ind) = ((ZCoords(iprot,tmpi,1)-ZCoords(iprot,tmpj,1))**2 
     .              + (ZCoords(iprot,tmpi,2)-ZCoords(iprot,tmpj,2))**2
     .              + (ZCoords(iprot,tmpi,3)-ZCoords(iprot,tmpj,3))**2)
               tmpj = tmpj + 1
               ind=ind+1
            enddo
            tmpi = tmpi + 1
            tmpj = zstrt2
         enddo

c     test cutoff criteria here before we do expensive operations
c     invrs currently contains the squares of the interatomic distances
         mindist = invrs(1)
         if ((mindist.gt.feather2).and.(cuttype.eq.RESIDUEMOLECULE)) then
c     as a heuristic check if invrs(1) is greater than the feather distance
c     only need to check other distances if true
            do ind=2,ui*uj
               if (invrs(ind) < mindist) mindist = invrs(ind)
               if (mindist.lt.feather2) goto 10
c     if we ever see a distance less than the feather distance then no need
c     to check more               
            enddo

            if ((mindist.gt.cut2).and.(cuttype.eq.RESIDUEMOLECULE)) then
               return
            endif
         endif
      
c     now transform invrs to actually be inverse interatomic separations
c     this is an expensive operation         
 10      do ind=1,ui*uj
            invrs(ind) = ONE / dsqrt ( invrs(ind) )
            invrs6(ind) = invrs(ind)**6
         enddo

c     Calculate energies
         ind = 1
         do i=zstrt1,zend1
            do j=zstrt2,zend2
               ipar1 = ZPar(iprot,i)
               ipar2 = ZPar(iprot,j)
               par = CLJPair(ipar1,ipar2)
c     Coulomb energy
               cnrg = cnrg + invrs(ind) * ChgIJ(par)
c     LJ energy               
               ljnrg = ljnrg + invrs6(ind)*(AIJ(par)*invrs6(ind) - BIJ(par))
               ind = ind + 1
            enddo
         enddo

         cnrgf = cnrg
         cnrgb = cnrg

         ljnrgf = ljnrg
         ljnrgb = ljnrg

      endif
      
      sepResResEnergy = cnrg + ljnrg
      
      return
      
      end
