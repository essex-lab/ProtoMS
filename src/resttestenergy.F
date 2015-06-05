      double precision function restTestEnergy()
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
      include 'quickpot.inc'
c###########################################################
c
c     This function returns the energy needs
c     for the REST Metropolis test, i.e.
c     all energy components involving solutes, scaled
c     differently than at a move 
c
c     (C) Samuel Genheden 2015
c
c###########################################################

      integer i,j
      double precision nrg
      double precision cnrg,ljnrg,mindist
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb
      double precision bndnrg,angnrg,dihnrg,ubynrg,nbcnrg,nbljnrg
      double precision bndnrgf,angnrgf,dihnrgf,ubynrgf,nbcnrgf,nbljnrgf
      double precision bndnrgb,angnrgb,dihnrgb,ubynrgb,nbcnrgb,nbljnrgb
      double precision scl,sclf,sclb
      double precision getLowScale,getHighScale
      double precision soluteIntraEnergy
      double precision soluteSoluteEnergy,soluteGCSEnergy
      double precision proteinSoluteEnergy,GCSoluteBulkEnergy
      double precision proteinGCSoluteEnergy
      double precision GCsoluteSolventEnergy
      double precision soluteSolventEnergy

      restTestEnergy = ZERO
      
c     Set scaling for test energy
      do i=1,7
        if (restlevel(i).eq.1) then
          restscaling(i) = 1.0d0
        else if (restlevel(i).eq.2) then
          restscaling(i) = 0.5d0
        else
          restscaling(i) = 0.0d0
        endif
      enddo
      
c     now do the same for all of the solutes...
      if (NSolutes.gt.0) then
        do i=1,NSolutes
          if (.not.IgnoreSolute(i)) then
            nrg = soluteIntraEnergy(i,
     .                bndnrg, angnrg, dihnrg, ubynrg,nbcnrg, nbljnrg,
     .                bndnrgb,angnrgb,dihnrgb,ubynrgb,nbcnrgb,nbljnrgb,
     .                bndnrgf,angnrgf,dihnrgf,ubynrgf,nbcnrgf,nbljnrgf)
            restTestEnergy = restTestEnergy + nrg
          endif
        enddo

        if (NSolutes.gt.1) then
          do i=1,NSolutes-1
            if (.not.IgnoreSolute(i)) then
              do j=i+1,NSolutes
                if (.not.IgnoreSolute(j)) then
c                 the scaling of solute-solute energies to account for dual topology
c                 solutes is done within the soluteSoluteEnergy function itself.
c                 This is because of the complexity of possibly having multiple dual topology
c                 solutes!
                  nrg = soluteSoluteEnergy(i,j,
     .                         cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist,
     .                         scl,sclf,sclb)
                  restTestEnergy = restTestEnergy + nrg                  
                endif
              enddo
            endif
          enddo
        endif

c       now the solute-protein energies...
        if (NProteins.gt.0) then
          do j=1,NSolutes
c           see if this is a dual topology solute
            if (DualTopologySolute(j).gt.0) then
              scl = getLowScale(Lambda)
              sclb = getLowScale(LambdaB)
              sclf = getLowScale(LambdaF)
            else if (DualTopologySolute(j).eq.0) then
              scl = ONE
              sclb = ONE
              sclf = ONE
            else
              scl = getHighScale(Lambda)
              sclb = getHighScale(LambdaB)
              sclf = getHighScale(LambdaF)
            endif
c          SG do theta scaling for fragments
           scl = scl*SoluteTheta(j) ! SG
           sclf = sclf*SoluteTheta(j) ! SG
           sclb = sclb*SoluteTheta(j) ! SG  
           
c          August 08. JM bugfix to remove harmless warnings when doing splitGBSA simulations with a dual
c          topology sims. The fwd/bwd energies should be the same, but for this to be true, the scale factors
c          must be the same at
           if (QuickMove) then
             sclb = scl
             sclf = scl
           endif
            if (.not.IgnoreSolute(j)) then
              do i=1,NProteins
                nrg = proteinSoluteEnergy(DONTCARE,i,j,
     .                      cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist,
     .                      scl,sclf,sclb)
               
c               scale the energies 
                nrg = nrg*scl
                cnrg = cnrg*scl
                ljnrg = ljnrg*scl
c               If the soft-core is on, the coulombic is scaled further here
                if (Solutesoftcore(j)) then
                  cnrg = cnrg * (scl)**CoulSoftPower
                endif                
c                totalEnergy = totalEnergy + nrg
                restTestEnergy = restTestEnergy + cnrg + ljnrg   
              enddo
            endif
          enddo
        endif
      endif

c########################################### 
c     now do the same for the GCS
      if (NGCSolutes.gt.0) then
c       calculate the bulk GCS energy
        nrg = GCSoluteBulkEnergy(
     .         cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
        restTestEnergy = restTestEnergy + nrg
      
c       now the solvent-GCS energy
        if (NSolvents.gt.0) then
c        add the total GCS-protein energy for protein i
          do i=1,NGCSolutes
            nrg = GCsoluteSolventEnergy(i,
     .              cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
            restTestEnergy = restTestEnergy + nrg         
          enddo
        endif

c       now the protein-GCS energy
        if (NProteins.gt.0) then
          do i=1,NProteins
c           add the total GCS-protein energy for protein i
            nrg = proteinGCSoluteEnergy(i,
     .              cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
            restTestEnergy = restTestEnergy + nrg
          enddo
        endif
        
c       now the solute-solvent energy
        if (NSolutes.gt.0) then
          do i=1,NSolutes
            if (.not.IgnoreSolute(i)) then

c             see whether we need to scale the energies to account for 
c             this being a dual topology solute
              if (DualTopologySolute(i).gt.0) then
                scl = getLowScale(Lambda)
                sclb = getLowScale(LambdaB)
                sclf = getLowScale(LambdaF)
              else if (DualTopologySolute(i).eq.0) then
                scl = ONE
                sclb = ONE
                sclf = ONE
              else
                scl = getHighScale(Lambda)
                sclb = getHighScale(LambdaB)
                sclf = getHighScale(LambdaF)
              endif
c             SG theta scaling for fragments
              scl = scl*SoluteTheta(i) ! SG
              sclf = sclf*SoluteTheta(i) ! SG
              sclb = sclb*SoluteTheta(i) ! SG  
              
              nrg = soluteGCSEnergy(DONTCARE,i,
     .              cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist,
     .              scl,sclf,sclb)

c             scale the energies
              nrg = nrg*scl
              cnrg = cnrg*scl
              ljnrg = ljnrg*scl
c             If the soft-core is on, the coulombic is scaled further here
              if (Solutesoftcore(i)) then
                cnrg = cnrg * (scl)**CoulSoftPower
              endif
              restTestEnergy = restTestEnergy + cnrg + ljnrg     
            endif
          enddo
        endif
      endif

c     now do the same for the solvent
      if (NSolvents.gt.0) then
    
c       now the solute-solvent energy
        if (NSolutes.gt.0) then
          do i=1,NSolutes
            if (.not.IgnoreSolute(i)) then

c             see whether we need to scale the energies to account for 
c             this being a dual topology solute
              if (DualTopologySolute(i).gt.0) then
                scl = getLowScale(Lambda)
                sclb = getLowScale(LambdaB)
                sclf = getLowScale(LambdaF)
              else if (DualTopologySolute(i).eq.0) then
                scl = ONE
                sclb = ONE
                sclf = ONE
              else
                scl = getHighScale(Lambda)
                sclb = getHighScale(LambdaB)
                sclf = getHighScale(LambdaF)
              endif
c             SG theta scaling for fragments
              scl = scl*SoluteTheta(i) ! SG
              sclf = sclf*SoluteTheta(i) ! SG
              sclb = sclb*SoluteTheta(i) ! SG  
              
              nrg = soluteSolventEnergy(DONTCARE,i,
     .              cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist,
     .              scl,sclf,sclb)
c             scale the energies
              nrg = nrg*scl
              cnrg = cnrg*scl
              ljnrg = ljnrg*scl
c             If the soft-core is on, the coulombic is scaled further here
              if (Solutesoftcore(i)) then
                cnrg = cnrg * (scl)**CoulSoftPower
              endif
              restTestEnergy = restTestEnergy + cnrg + ljnrg     
            endif
          enddo
        endif
      endif

      do i=1,7
        if (restlevel(i).eq.1) then
          restscaling(i) = restfac1
        else if (restlevel(i).eq.2) then
          restscaling(i) = restfac2
        else
          restscaling(i) = 1.0d0
        endif
      enddo

      return

      end
