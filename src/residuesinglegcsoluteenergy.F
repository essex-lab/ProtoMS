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
      integer i,j,ui,uj,tmpi,tmpj,ind,par,nat,resstrt,resend
      double precision sclfac
      double precision mindist
      include 'nbenergyvar2.inc'
      double precision charge,sigma,epslon      
      double precision invr2, tmp2, tmp3, tmp
c     Size of these arrays imply an implicit limit on size of gcsolute
c     which doesn't seem to be consistently handled. The practical
c     limit is determined by the definition of ZGCSPar in solvents.inc at
c     50 atoms.
      double precision invrs((MAXSCATOMS+3)*50),invrs6((MAXSCATOMS+3)*50)

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
      
      if (ddd.or.softcore.or.(boundary.eq.PERIODICBOUNDARIES)) then
c     run unoptmised version of nbenergy as this can
c     handle softcores, implicit solvent and periodic boundary conditions
         do i=resstrt,resend
            funcat1 = i
            
            do j=1,nat
               funcat2 = j

               include 'nbinitial.inc'
               include 'nbprogcs.inc'
               include 'nbenergy3.inc'
               
               nrg = nrg + funcnbnrg
               cnrg = cnrg + funccnrg
c     PRINT*, "cnrg" ,cnrg
               ljnrg = ljnrg + funcljnrg
c     PRINT*, "ljnrg" ,ljnrg
               cnrgf = cnrgf + funccnrgf
               cnrgb = cnrgb + funccnrgb
               ljnrgf = ljnrgf + funcljnrgf
               ljnrgb = ljnrgb + funcljnrgb
               if (funcmindist.lt.mindist) mindist = funcmindist
            enddo
         enddo
      else
c     optimised version of nbenergy cannot handle
c     softcores, implicit solvent or lambda scaling

         thet = thetaparam (igcs)
         ui=resend-resstrt+1
         uj=nat

         tmpi=resstrt
         tmpj=1
c        PERIODIC BOUNDARIES!!!!!!!!
         ind=1
         do i=1,ui
            do j=1,uj
               invrs(ind) = ((ZCoords(iprot,tmpi,1)-ZGCSCoords(igcs,tmpj,1))**2
     .              + (ZCoords(iprot,tmpi,2)-ZGCSCoords(igcs,tmpj,2))**2
     .              + (ZCoords(iprot,tmpi,3)-ZGCSCoords(igcs,tmpj,3))**2)
               tmpj = tmpj + 1
               ind=ind+1
            enddo
            tmpi = tmpi + 1
            tmpj = 1
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
c$$$  goto 20
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
         do i=resstrt,resend
            do j=1,nat
               ipar1 = ZPar(iprot,i)
               ipar2 = ZGCSPar(igcs,j)
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

         nrg = cnrg + ljnrg
      endif
      
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

      if (doREST) then
        cnrg = cnrg * restscaling(6)
        ljnrg = ljnrg * restscaling(6)
        nrg = cnrg + ljnrg
        cnrgf = cnrgf * restscaling(6)
        ljnrgf = ljnrgf * restscaling(6)
        cnrgb = cnrgb * restscaling(6)
        ljnrgb = ljnrgb * restscaling(6)
      endif

      residueSingleGCsoluteEnergy = nrg
c      PRINT*,"residue-gc nrg" ,residueSingleGCsoluteEnergy
      return
      end
