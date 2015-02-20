      double precision function pmeEnergy(selfnrg,exclnrg,meshnrg,chnrg,selfsvnnrg,exclsvnnrg,initialize)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'simparams.inc'
      include 'solvents.inc'
      include 'movelist.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'energies.inc'
c###########################################################
c
c     This function calculates the self, and
c     exclusion energy that is necessary for a particle-
c     mesh Ewald calculations, it also call an external
c     C++ routine to calculate the repciprocal/mesh energy
c
c     (C) Samuel Genheden        27 September 2013
c
c###########################################################

      integer MAXPMEATOMS 
      parameter (MAXPMEATOMS = 50000)
      integer natoms0,natoms
      double precision x(MAXPMEATOMS),y(MAXPMEATOMS),z(MAXPMEATOMS),charges(MAXPMEATOMS)
      double precision cutoff_nm,box_nm(3),pme_alpha_nm 

      double precision selfnrg, exclnrg,meshnrg,chnrg
      double precision selfnrg1,exclnrg1,meshnrg1,chnrg1
      double precision selfnrg2,exclnrg2,meshnrg2,chnrg2
      double precision selfsvnnrg,exclsvnnrg
      double precision chsum,chsum1,chsum2,fac,vol

      integer i,j,ndual
      double precisionlam1,lam2
      double precision mesh_energy,getLowScale,getHighScale
      logical initialize

c     Convert cutoff and box dimensions to nm
      cutoff_nm = cut/10.0d0
      box_nm(1) = BoxDimension(1)/10.0d0
      box_nm(2) = BoxDimension(2)/10.0d0
      box_nm(3) = BoxDimension(3)/10.0d0
      pme_alpha_nm = pme_alpha*10.0d0

c     Pre-factor for charge-correction
      vol = BoxDimension(1)*BoxDimension(2)*BoxDimension(3)
      fac = FAC4PIEPS0*PI/(2.0d0*vol*vol*pme_alpha*pme_alpha)
      chsum = ZERO
      chsum1 = ZERO
      chsum2 = ZERO

c     No atoms from beginning
      natoms = 0

c     Zero all energies
      selfnrg1 = ZERO
      selfnrg2 = ZERO
      selfnrg  = ZERO
      exclnrg1 = ZERO
      exclnrg2 = ZERO
      exclnrg  = ZERO
      meshnrg1 = ZERO
      meshnrg2 = ZERO
      meshnrg  = ZERO

c     Calculate contribution from the solvent molecules      
      if (initialize) then
        selfsvnnrg = ZERO
        exclsvnnrg = ZERO 
      endif      
      call solventPMEEnergy(exclsvnnrg,selfsvnnrg,charges,x,y,z,natoms,MAXPMEATOMS,initialize)
      selfnrg1 = selfsvnnrg
      exclnrg1 = exclsvnnrg

c     Calculate contribution from all the proteins
      call proteinPMEEnergy(exclnrg1,selfnrg1,charges,x,y,z,natoms,MAXPMEATOMS)          

c     Calculate contribution from all solutes that do not have a dual-topology partner
      ndual = 0
      do i=1,NSolutes
        if (DualTopologySolute(i).eq.0) then
          call solutePMEEnergy(i,exclnrg1,selfnrg1,charges,x,y,z,natoms,MAXPMEATOMS)
        else
          ndual = ndual + 1
        endif
      enddo           

c     Compute partial charge-sum
      do i=1,natoms
        chsum = chsum + charges(i)
      enddo

c     Simple do only one mesh calculations
      if (ndual.eq.0) then
        chnrg = -vol*chsum*chsum*fac
        meshnrg1 = mesh_energy(natoms,box_nm,x,y,z,charges,cutoff_nm,pme_tol,pme_alpha_nm,0)
        meshnrg1 = meshnrg1/4.184d0
        selfnrg = selfnrg1
        exclnrg = exclnrg1
        meshnrg = meshnrg1
c     Have to deal with dual topology solutes
      else
c       This stores the number of atoms loaded before dual-topology solutes
        natoms0 = natoms

c       Set energy of state2 (V1) to energy of state1 (V0)
        selfnrg2 = selfnrg1
        exclnrg2 = exclnrg1

c       Charge sums
        chsum1 = chsum
        chsum2 = chsum

c       Load one solute, V0
        do i=1,NSolutes
          if (DualTopologySolute(i).gt.0)
     .       call solutePMEEnergy(i,exclnrg1,selfnrg1,charges,x,y,z,natoms)
        enddo        
        do i=natoms0+1,natoms
          chsum1 = chsum1 + charges(i)
        enddo

c        write(6,*)natoms,exclnrg1+selfnrg1+meshnrg1,exclnrg1,selfnrg1,meshnrg1

        meshnrg1 = mesh_energy(natoms,box_nm,x,y,z,charges,cutoff_nm,pme_tol,pme_alpha_nm,0)
        meshnrg1 = meshnrg1/4.184d0
           
c        write(6,*)natoms,exclnrg1+selfnrg1+meshnrg1,exclnrg1,selfnrg1,meshnrg1

c       Load the other solute, V1
c       This will overwrite charges and coordinates for solutes in V0
        natoms = natoms0
        do i=1,NSolutes
          if (DualTopologySolute(i).lt.0) 
     .       call solutePMEEnergy(i,exclnrg2,selfnrg2,charges,x,y,z,natoms)
        enddo        
        do i=natoms0+1,natoms
          chsum2 = chsum2 + charges(i)
        enddo

        meshnrg2 = mesh_energy(natoms,box_nm,x,y,z,charges,cutoff_nm,pme_tol,pme_alpha_nm,0)
        meshnrg2 = meshnrg2/4.184d0

c        write(6,*)natoms,exclnrg2+selfnrg2+meshnrg2,exclnrg2,selfnrg2,meshnrg2

        lam1 = getLowScale(Lambda)
        lam2 = getHighScale(Lambda)
        selfnrg = lam1*selfnrg1+lam2*selfnrg2
        exclnrg = lam1*exclnrg1+lam2*exclnrg2
        meshnrg = lam1*meshnrg1+lam2*meshnrg2
        chnrg   = -vol*(lam1*chsum1+lam2*chsum2)*fac

      endif

      pmeEnergy = selfnrg + exclnrg + meshnrg + chnrg

      return

      end function
