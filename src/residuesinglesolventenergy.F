      double precision function residueSingleSolventEnergy(iprot,ires,isvn,
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
c####################################################################
c
c     This calculates the interaction energy between residue ires
c     of protein iprot, and solvent molecule isvn
c
c     (C) Christopher Woods    29 January 2002
c
c     Extensive rewrite to improve optimisation of code
c
c     (C) Chris Cave-Ayland March 2015     
c     
c####################################################################

      integer iprot,ires,isvn
      double precision cnrg,ljnrg,nrg
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb
      double precision distx,disty,distz,dist2
      integer i,j,nat,resstrt,resend,istart
      double precision sclfac
      double precision mindist
      include 'nbenergyvar2.inc'
      double precision charge,sigma,epslon      
      double precision getSolventRadius
c     temporary cheat
      double precision scl,sclf,sclb
      double precision invrs(120), invrs6(120)
      integer ind, par, tmpj, ui, uj
      
      scl = 1.0d0
      sclf = 1.0d0
      sclb = 1.0d0
c     end of cheat
      residueSingleSolventEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      nrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      mindist = 99900.0
      nat = NSvnAtoms(isvn)
      softcore = .false.
      ddd = .false.
      if (nat.le.0) return

c     return if either the residue or solvent has been done before
      if (DoneResNrg(iprot,ires) .or. DoneSvnNrg(isvn)) return

c     now check the distance between the residue centre and centre
c     of the solvent molecule - cut if it greater than cutoff+residueRadius
      distx = dabs(ResCentre(iprot,ires,1)-SvnCentre(isvn,1))
      disty = dabs(ResCentre(iprot,ires,2)-SvnCentre(isvn,2))
      distz = dabs(ResCentre(iprot,ires,3)-SvnCentre(isvn,3))
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
        if (dist2 .gt. (cut+ResRadius(iprot,ires))**2) return
      else
        if (dist2 .gt. (cut+ResRadius(iprot,ires)+getSolventRadius(isvn))**2) return
      endif     
 
c     assume residue runs concurrently from N to O
      resstrt = BBZid(iprot,ires,1)
      resend = BBZid(iprot,ires,4)
      funciprot1 = iprot
      funcisvn2 = isvn
      if (QuickMove .and. DDDon) ddd = .true.
      if (ddd) then
c     use unoptimised nbenergy routines that can handle implicit solvent
         do i=resstrt,resend
            funcat1 = i
            if ((SvnOptimization(isvn).eq.OPTTIP3P.or.SvnOptimization(isvn).eq.OPTTIP4P)) then
            
c     Optimized for TIP4P-TIP4P interactions
            
c     LJ on oxygens only
              funcat2 = 1
              include 'nbinitial.inc'
              include 'nbprosvn.inc'
              include 'nbenergy_lj.inc'         
              
              nrg = nrg + funcnbnrg
              ljnrg = ljnrg + funcljnrg
              ljnrgf = ljnrgf + funcljnrgf
              ljnrgb = ljnrgb + funcljnrgb
              if (funcmindist.lt.mindist) mindist = funcmindist
              
c     Charges on the rest of the particles, i.e., between hydrogens and extra point particle
              if (SvnOptimization(isvn).eq.OPTTIP3P) then
                 istart = 1
              else
                 istart = 2
              endif
              do j=istart,nat
                 funcat2 = j
                 
                 include 'nbinitial.inc'
                 include 'nbprosvn.inc'
                 include 'nbenergy_coul.inc'
                 
                 nrg = nrg + funcnbnrg
                 cnrg = cnrg + funccnrg
                 cnrgf = cnrgf + funccnrgf
                 cnrgb = cnrgb + funccnrgb
                 
                 if (funcmindist.lt.mindist) mindist = funcmindist
              enddo
           else              
              do j=1,nat
                 funcat2 = j
                 
                 include 'nbinitial.inc'
                 include 'nbprosvn.inc'
                 include 'nbenergy2.inc'
                 
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
      else
c     use optimised nbenergy routines

         ui = nat
         uj = resend-resstrt+1
         
         tmpj = resstrt
         
         ind=1
         do i=1,ui
            do j=1,uj
               invrs(ind) = ((ZCoords(iprot,tmpj,1)-ZSvnCoords(isvn,i,1))**2 
     .              + (ZCoords(iprot,tmpj,2)-ZSvnCoords(isvn,i,2))**2
     .              + (ZCoords(iprot,tmpj,3)-ZSvnCoords(isvn,i,3))**2)
               tmpj = tmpj + 1
               ind=ind+1
            enddo
            tmpj = resstrt
         enddo
         
c     test cutoff criteria here before we do expensive operations
c     invrs currently contains the squares of the interatomic distances
         mindist = invrs(1)
         if ((mindist.gt.feather2).and.(cuttype.eq.RESIDUEMOLECULE)) then
c     as a heuristic, check if invrs(1) is greater than the feather distance
c     only need to check other distances if true
            do ind=2,ui*uj
               if (invrs(ind) < mindist) mindist = invrs(ind)
               if (mindist.lt.feather2) goto 10
c     if we ever see a distance less than the feather distance then no need
c     to check more               
            enddo
            if ((mindist.gt.cut2).and.(cuttype.eq.RESIDUEMOLECULE)) then
               nrg = ZERO
               cnrg = ZERO
               ljnrg = ZERO
               cnrgf = ZERO
               cnrgb = ZERO
               ljnrgf = ZERO
               ljnrgb = ZERO
               sclfac = ZERO
               ResidueSingleSolventEnergy = ZERO
               return
            endif
         endif

c     now transform invrs to actually be inverse interatomic separations
c     this is an expensive operation         
 10      do ind=1,ui*uj
            invrs(ind) = ONE / dsqrt ( invrs(ind) )
         enddo
         
c     we can optimize here as we only need ivrs6s for O if TIP3P / TIP4P
         if ((SvnOptimization(isvn).eq.OPTTIP3P).or.(SvnOptimization(isvn).eq.OPTTIP4P)) then
            do ind=1,uj
               invrs6(ind) = invrs(ind)**6
            enddo
         else
            do ind=1,ui*uj
               invrs6(ind) = invrs(ind)**6
            enddo
         endif
         

         if ((SvnOptimiZation(isvn).eq.OPTTIP3P).or.(SvnOptimiZation(isvn).eq.OPTTIP4P)) then
            if (SvnOptimiZation(isvn).eq.OPTTIP3P) then
c     for tip3p we do both lj and colomb for the oxygen (at1)
               ind = 1
               do j=resstrt,resend
                  ipar1 = ZPar(iprot,j)
                  ipar2 = ZSvnPar(isvn,1)
                  par = CLJPair(ipar1,ipar2)
c     Coulomb energy
                  cnrg = cnrg + invrs(ind) * ChgIJ(par)
c     LJ energy               
                  ljnrg = ljnrg + invrs6(ind)*(AIJ(par)*invrs6(ind) - BIJ(par))
                  ind = ind + 1
               enddo

            elseif (SvnOptimiZation(isvn).eq.OPTTIP4P) then
c     for tip4p we do lj for the oxygen (at1)
               ind = 1
               do j=resstrt,resend
                  ipar1 = ZPar(iprot,j)
                  ipar2 = ZSvnPar(isvn,1)
                  par = CLJPair(ipar1,ipar2)
c     LJ energy               
                  ljnrg = ljnrg + invrs6(ind)*(AIJ(par)*invrs6(ind) - BIJ(par))
                  ind = ind + 1
               enddo
            endif               

c     now do Coulomb energy for the rest of the atoms
            ind = uj + 1
            do i=2,nat
               do j=resstrt,resend
                  ipar1 = ZPar(iprot,j)
                  ipar2 = ZSvnPar(isvn,i)
                  par = CLJPair(ipar1,ipar2)

                  cnrg = cnrg + invrs(ind) * ChgIJ(par)
                  ind = ind + 1
               enddo
            enddo

         else
c     anything else and we do lj and coulomb for all
            ind = 1
            do i=1,nat
               do j=resstrt,resend
                  ipar1 = ZPar(iprot,j)
                  ipar2 = ZSvnPar(isvn,i)
                  par = CLJPair(ipar1,ipar2)
c     Coulomb energy
                  cnrg = cnrg + invrs(ind) * ChgIJ(par)
c     LJ energy               
                  ljnrg = ljnrg + invrs6(ind)*(AIJ(par)*invrs6(ind) - BIJ(par))
                  ind = ind + 1
               enddo
            enddo

         endif
         cnrgf = cnrg
         cnrgb = cnrg

         ljnrgf = ljnrg
         ljnrgb = ljnrg
         nrg = ljnrg + cnrg
         
      endif 
      if (cuttype.eq.RESIDUEMOLECULE) then
c       doing protein-residue - solvent molecule cutting
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
      
      residueSingleSolventEnergy = nrg
  
      return
      end
