      double precision function totalEnergy(savenrg,nrgb,nrgf)
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
c     This function returns the total energy of the 
c     entire system! Obviously this is quite an expensive
c     thing to do, though it should rarely need to be done...
c
c     If savenrg=OLD then the components etc. are saved as the
c     'old' energy. If savenrg=NEW they are saved as the 'new'
c     energy. Otherwise the components are lost to the ether!
c
c
c     (C) Christopher Woods    7 February 2002
c
c     Modified by Samuel Genheden 2012-2014
c
c###########################################################

      integer savenrg
      integer i,j
      double precision nrg,nrgb,nrgf
      double precision cnrg,ljnrg,mindist
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb
      double precision bndnrg,angnrg,dihnrg,ubynrg,nbcnrg,nbljnrg
      double precision bndnrgf,angnrgf,dihnrgf,ubynrgf,nbcnrgf,nbljnrgf
      double precision bndnrgb,angnrgb,dihnrgb,ubynrgb,nbcnrgb,nbljnrgb
      double precision scl,sclf,sclb
      double precision getLowScale,getHighScale
      double precision proteinIntraEnergy
      double precision proteinProteinEnergy
      double precision soluteIntraEnergy
      double precision soluteSoluteEnergy,soluteGCSEnergy
      double precision proteinSoluteEnergy,GCSoluteBulkEnergy
      double precision solventSolventBulkEnergy,proteinGCSoluteEnergy
      double precision proteinSolventenergy,GCsoluteSolventEnergy
      double precision soluteSolventEnergy

      totalEnergy = ZERO
      nrgb = ZERO
      nrgf = ZERO

      if (NProteins.gt.0) then
c       loop over all proteins and get them to return thier
c       intramolecular energy...
        do i=1,NProteins
          nrg = proteinIntraEnergy(i,
     .            bndnrg, angnrg, dihnrg, ubynrg, nbcnrg, nbljnrg,
     .            bndnrgb,angnrgb,dihnrgb,ubynrgb,nbcnrgb,nbljnrgb,
     .            bndnrgf,angnrgf,dihnrgf,ubynrgf,nbcnrgf,nbljnrgf)
          totalEnergy = totalEnergy + nrg
          nrgf = nrgf + bndnrgf+angnrgf+dihnrgf+ubynrgf+nbcnrgf+nbljnrgf
          nrgb = nrgb + bndnrgb+angnrgb+dihnrgb+ubynrgb+nbcnrgb+nbljnrgb

          if (savenrg.eq.OLD) then
            OldProBnd(i) = bndnrg
            OldProAng(i) = angnrg
            OldProDih(i) = dihnrg
            OldProUBy(i) = ubynrg
            OldProNBC(i) = nbcnrg
            OldProNBLJ(i) = nbljnrg
            OldProBndF(i) = bndnrgf
            OldProAngF(i) = angnrgf
            OldProDihF(i) = dihnrgf
            OldProUByF(i) = ubynrgf
            OldProNBCF(i) = nbcnrgf
            OldProNBLJF(i) = nbljnrgf
            OldProBndB(i) = bndnrgb
            OldProAngB(i) = angnrgb
            OldProDihB(i) = dihnrgb
            OldProUByB(i) = ubynrgb
            OldProNBCB(i) = nbcnrgb
            OldProNBLJB(i) = nbljnrgb
          else if (savenrg.eq.NEW) then
            NewProBnd(i) = bndnrg
            NewProAng(i) = angnrg
            NewProDih(i) = dihnrg
            NewProUBy(i) = ubynrg
            NewProNBC(i) = nbcnrg
            NewProNBLJ(i) = nbljnrg
            NewProBndF(i) = bndnrgf
            NewProAngF(i) = angnrgf
            NewProDihF(i) = dihnrgf
            NewProUByF(i) = ubynrgf
            NewProNBCF(i) = nbcnrgf
            NewProNBLJF(i) = nbljnrgf
            NewProBndB(i) = bndnrgb
            NewProAngB(i) = angnrgb
            NewProDihB(i) = dihnrgb
            NewProUByB(i) = ubynrgb
            NewProNBCB(i) = nbcnrgb
            NewProNBLJB(i) = nbljnrgb
          endif        
        enddo

      
c       now loop over all pairs of proteins and calculate thier
c       protein-protein energy...
        if (NProteins.gt.1) then
          do i=1,NProteins-1
            do j=i+1,NProteins
              totalEnergy = totalEnergy + proteinProteinEnergy(i,j,
     .                   cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
              nrgf = nrgf + cnrgf + ljnrgf
              nrgb = nrgb + cnrgb + ljnrgb
              if (savenrg.eq.OLD) then
                OldProProC(i,j) = cnrg
                OldProProLJ(i,j) = ljnrg
                OldProProC(j,i) = cnrg
                OldProProLJ(j,i) = ljnrg
                OldProProCF(i,j) = cnrgf
                OldProProLJF(i,j) = ljnrgf
                OldProProCF(j,i) = cnrgf
                OldProProLJF(j,i) = ljnrgf
                OldProProCB(i,j) = cnrgb
                OldProProLJB(i,j) = ljnrgb
                OldProProCB(j,i) = cnrgb
                OldProProLJB(j,i) = ljnrgb
                OldProProDist(j,i) = mindist
                OldProProDist(i,j) = mindist
              else if (savenrg.eq.NEW) then
                NewProProC(i,j) = cnrg
                NewProProLJ(i,j) = ljnrg
                NewProProC(j,i) = cnrg
                NewProProLJ(j,i) = ljnrg
                NewProProCF(i,j) = cnrgf
                NewProProLJF(i,j) = ljnrgf
                NewProProCF(j,i) = cnrgf
                NewProProLJF(j,i) = ljnrgf
                NewProProCB(i,j) = cnrgb
                NewProProLJB(i,j) = ljnrgb
                NewProProCB(j,i) = cnrgb
                NewProProLJB(j,i) = ljnrgb
                NewProProDist(j,i) = mindist
                NewProProDist(i,j) = mindist
              endif
            enddo
          enddo
        endif
      endif
      
c     now do the same for all of the solutes...
      if (NSolutes.gt.0) then
        do i=1,NSolutes
          if (.not.IgnoreSolute(i)) then
            nrg = soluteIntraEnergy(i,
     .                bndnrg, angnrg, dihnrg, ubynrg,nbcnrg, nbljnrg,
     .                bndnrgb,angnrgb,dihnrgb,ubynrgb,nbcnrgb,nbljnrgb,
     .                bndnrgf,angnrgf,dihnrgf,ubynrgf,nbcnrgf,nbljnrgf)

c           see if the solute is part of a dual topology solute
            if (DualTopologySolute(i).gt.0) then
c             this is the low dual topology solute
              scl = getLowScale(Lambda)
              sclb = getLowScale(LambdaB)
              sclf = getLowScale(LambdaF)
            else if (DualTopologySolute(i).eq.0) then
c             this is not a dual topology solute
              scl = ONE
              sclf = ONE
              sclb = ONE
            else
c             this is the high dual topology solute
              scl = getHighScale(Lambda)
              sclb = getHighScale(LambdaB)
              sclf = getHighScale(LambdaF)
            endif

c           The intramolecular energy is not scaled
            nrg = nrg
            bndnrg = bndnrg
            angnrg = angnrg
            dihnrg = dihnrg
            ubynrg = ubynrg
            nbcnrg = nbcnrg
            nbljnrg = nbljnrg
            bndnrgb = bndnrgb
            angnrgb = angnrgb
            dihnrgb = dihnrgb
            ubynrgb = ubynrgb
            nbcnrgb = nbcnrgb
            nbljnrgb = nbljnrgb
            bndnrgf = bndnrgf
            angnrgf = angnrgf
            dihnrgf = dihnrgf
            ubynrgf = ubynrgf
            nbcnrgf = nbcnrgf
            nbljnrgf = nbljnrgf
c           No softcore is applied to the intramolecular energy
            totalEnergy = totalEnergy + nrg
            nrgf = nrgf + bndnrgf+angnrgf+dihnrgf+ubynrgf+nbcnrgf+nbljnrgf
            nrgb = nrgb + bndnrgb+angnrgb+dihnrgb+ubynrgb+nbcnrgb+nbljnrgb

            if (savenrg.eq.OLD) then
              OldSolBnd(i) = bndnrg
              OldSolAng(i) = angnrg
              OldSolDih(i) = dihnrg
              OldSolUBy(i) = ubynrg
              OldSolNBC(i) = nbcnrg
              OldSolNBLJ(i) = nbljnrg
              OldSolBndF(i) = bndnrgf
              OldSolAngF(i) = angnrgf
              OldSolDihF(i) = dihnrgf
              OldSolUByF(i) = ubynrgf
              OldSolNBCF(i) = nbcnrgf
              OldSolNBLJF(i) = nbljnrgf
              OldSolBndB(i) = bndnrgb  
              OldSolAngB(i) = angnrgb
              OldSolDihB(i) = dihnrgb
              OldSolUByB(i) = ubynrgb
              OldSolNBCB(i) = nbcnrgb
              OldSolNBLJB(i) = nbljnrgb
            else if (savenrg.eq.NEW) then
              NewSolBnd(i) = bndnrg
              NewSolAng(i) = angnrg
              NewSolDih(i) = dihnrg
              NewSolUBy(i) = ubynrg
              NewSolNBC(i) = nbcnrg
              NewSolNBLJ(i) = nbljnrg
              NewSolBndF(i) = bndnrgf
              NewSolAngF(i) = angnrgf
              NewSolDihF(i) = dihnrgf
              NewSolUByF(i) = ubynrgf
              NewSolNBCF(i) = nbcnrgf
              NewSolNBLJF(i) = nbljnrgf
              NewSolBndB(i) = bndnrgb
              NewSolAngB(i) = angnrgb
              NewSolDihB(i) = dihnrgb
              NewSolUByB(i) = ubynrgb
              NewSolNBCB(i) = nbcnrgb
              NewSolNBLJB(i) = nbljnrgb
            endif        
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
                  totalEnergy = totalEnergy + nrg
                  nrgf = nrgf + cnrgf + ljnrgf
                  nrgb = nrgb + cnrgb + ljnrgb
                  if (savenrg.eq.OLD) then
                    OldSolSolC(j,i) = cnrg
                    OldSolSolLJ(j,i) = ljnrg
                    OldSolSolC(i,j) = cnrg
                    OldSolSolLJ(i,j) = ljnrg
                    OldSolSolCF(j,i) = cnrgf
                    OldSolSolLJF(j,i) = ljnrgf
                    OldSolSolCF(i,j) = cnrgf
                    OldSolSolLJF(i,j) = ljnrgf
                    OldSolSolCB(j,i) = cnrgb
                    OldSolSolLJB(j,i) = ljnrgb
                    OldSolSolCB(i,j) = cnrgb
                    OldSolSolLJB(i,j) = ljnrgb
                    OldSolSolDist(j,i) = mindist
                    OldSolSolDist(i,j) = mindist
                  else if (savenrg.eq.NEW) then
                    NewSolSolC(j,i) = cnrg
                    NewSolSolLJ(j,i) = ljnrg
                    NewSolSolC(i,j) = cnrg
                    NewSolSolLJ(i,j) = ljnrg
                    NewSolSolCF(j,i) = cnrgf
                    NewSolSolLJF(j,i) = ljnrgf
                    NewSolSolCF(i,j) = cnrgf
                    NewSolSolLJF(i,j) = ljnrgf
                    NewSolSolCB(j,i) = cnrgb
                    NewSolSolLJB(j,i) = ljnrgb
                    NewSolSolCB(i,j) = cnrgb
                    NewSolSolLJB(i,j) = ljnrg
                    NewSolSolDist(j,i) = mindist
                    NewSolSolDist(i,j) = mindist
                  endif                    
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
                cnrgf = cnrgf*sclf
                ljnrgf = ljnrgf*sclf
                cnrgb = cnrgb*sclb
                ljnrgb = ljnrgb*sclb
c               If the soft-core is on, the coulombic is scaled further here
                if (Solutesoftcore(j)) then
                  cnrg = cnrg * (scl)**CoulSoftPower
                  cnrgf = cnrgf * (sclf)**CoulSoftPower
                  cnrgb = cnrgb * (sclb)**CoulSoftPower
                endif                
                nrgf = nrgf + cnrgf + ljnrgf
                nrgb = nrgb + cnrgb + ljnrgb
c                totalEnergy = totalEnergy + nrg
                totalEnergy = totalEnergy + cnrg + ljnrg
                if (savenrg.eq.OLD) then
                  OldProSolC(i,j) = cnrg
                  OldProSolLJ(i,j) = ljnrg
                  OldProSolCF(i,j) = cnrgf
                  OldProSolLJF(i,j) = ljnrgf
                  OldProSolCB(i,j) = cnrgb
                  OldProSolLJB(i,j) = ljnrgb
                  OldProSolDist(i,j) = mindist
                else if (savenrg.eq.NEW) then
                  NewProSolC(i,j) = cnrg
                  NewProSolLJ(i,j) = ljnrg
                  NewProSolCF(i,j) = cnrgf
                  NewProSolLJF(i,j) = ljnrgf
                  NewProSolCB(i,j) = cnrgb
                  NewProSolLJB(i,j) = ljnrgb
                  NewProSolDist(i,j) = mindist
                endif      
              enddo
            endif
          enddo
        endif
      endif

c      write(6,*)totalEnergy,nrgf,nrgb


c########################################### 
c     now do the same for the GCS
      if (NGCSolutes.gt.0) then
c       calculate the bulk GCS energy
        nrg = GCSoluteBulkEnergy(
     .         cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
        nrgf = nrgf + cnrgf + ljnrgf
        nrgb = nrgb + cnrgb + ljnrgb
        totalEnergy = totalEnergy + nrg
        
        if (savenrg.eq.OLD) then
          OldGCSGCSC = cnrg
          OldGCSGCSLJ = ljnrg
          OldGCSGCSCF = cnrgf
          OldGCSGCSLJF = ljnrgf
          OldGCSGCSCB = cnrgb
          OldGCSGCSLJB = ljnrgb
        else if (savenrg.eq.NEW) then
          NewGCSGCSC = cnrg
          NewGCSGCSLJ = ljnrg
          NewGCSGCSCF = cnrgf
          NewGCSGCSLJF = ljnrgf
          NewGCSGCSCB = cnrgb
          NewGCSGCSLJB = ljnrgb
        endif
      
c       now the solvent-GCS energy
        if (NSolvents.gt.0) then
c        add the total GCS-protein energy for protein i
          do i=1,NGCSolutes
            nrg = GCsoluteSolventEnergy(i,
     .              cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
            nrgf = nrgf + cnrgf + ljnrgf
            nrgb = nrgb + cnrgb + ljnrgb
            totalEnergy = totalEnergy + nrg
            if (savenrg.eq.OLD) then
              OldGCSSvnC(i) = cnrg
              OldGCSSvnLJ(i) = ljnrg
              OldGCSSvnCF(i) = cnrgf
              OldGCSSvnLJF(i) = ljnrgf
              OldGCSSvnCB(i) = cnrgb
              OldGCSSvnLJB(i) = ljnrgb
            else if (savenrg.eq.NEW) then
              NewGCSSvnC(i) = cnrg
              NewGCSSvnLJ(i) = ljnrg
              NewGCSSvnCF(i) = cnrgf
              NewGCSSvnLJF(i) = ljnrgf
              NewGCSSvnCB(i) = cnrgb
              NewGCSSvnLJB(i) = ljnrgb
            endif         
          enddo
        endif

c       write(6,*)totalEnergy,nrgf,nrgb

c       now the protein-GCS energy
        if (NProteins.gt.0) then
          do i=1,NProteins
c           add the total GCS-protein energy for protein i
            nrg = proteinGCSoluteEnergy(i,
     .              cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
            nrgf = nrgf + cnrgf + ljnrgf
            nrgb = nrgb + cnrgb + ljnrgb
            totalEnergy = totalEnergy + nrg
            if (savenrg.eq.OLD) then
              OldProGCSC(i) = cnrg
              OldProGCSLJ(i) = ljnrg
              OldProGCSCF(i) = cnrgf
              OldProGCSLJF(i) = ljnrgf
              OldProGCSCB(i) = cnrgb
              OldProGCSLJB(i) = ljnrgb
            else if (savenrg.eq.NEW) then
              NewProGCSC(i) = cnrg
              NewProGCSLJ(i) = ljnrg
              NewProGCSCF(i) = cnrgf
              NewProGCSLJF(i) = ljnrgf
              NewProGCSCB(i) = cnrgb
              NewProGCSLJB(i) = ljnrgb
            endif
          enddo
        endif

c      write(6,*)totalEnergy,nrgf,nrgb
        
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

c      write(6,*)"H",nrg,cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf

c             scale the energies
              nrg = nrg*scl
              cnrg = cnrg*scl
              ljnrg = ljnrg*scl
              cnrgf = cnrgf*sclf
              ljnrgf = ljnrgf*sclf
              cnrgb = cnrgb*sclb
              ljnrgb = ljnrgb*sclb
c      write(6,*)ljnrg,ljnrgf,ljnrgb,nrg,cnrg+ljnrg
c             If the soft-core is on, the coulombic is scaled further here
              if (Solutesoftcore(i)) then
                cnrg = cnrg * (scl)**CoulSoftPower
                cnrgf = cnrgf * (sclf)**CoulSoftPower
                cnrgb = cnrgb * (sclb)**CoulSoftPower
              endif
              nrgf = nrgf + cnrgf + ljnrgf
              nrgb = nrgb + cnrgb + ljnrgb
c             totalEnergy = totalEnergy + nrg
              totalEnergy = totalEnergy + cnrg + ljnrg
              
              if (savenrg.eq.OLD) then
                OldSolGCSC(i) = cnrg
                OldSolGCSLJ(i) = ljnrg
                OldSolGCSCF(i) = cnrgf
                OldSolGCSLJF(i) = ljnrgf
                OldSolGCSCB(i) = cnrgb
                OldSolGCSLJB(i) = ljnrgb
              else if (savenrg.eq.NEW) then
                NewSolGCSC(i) = cnrg
                NewSolGCSLJ(i) = ljnrg
                NewSolGCSCF(i) = cnrgf
                NewSolGCSLJF(i) = ljnrgf
                NewSolGCSCB(i) = cnrgb
                NewSolGCSLJB(i) = ljnrgb
              endif      
            endif
          enddo
        endif
      endif
      
c####################################



c     now do the same for the solvent
      if (NSolvents.gt.0) then
c       calculate the bulk solvent-solvent energy
        nrg = solventSolventBulkEnergy(
     .         cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
        nrgf = nrgf + cnrgf + ljnrgf
        nrgb = nrgb + cnrgb + ljnrgb
        totalEnergy = totalEnergy + nrg
        
        if (savenrg.eq.OLD) then
          OldSvnSvnC = cnrg
          OldSvnSvnLJ = ljnrg
          OldSvnSvnCF = cnrgf
          OldSvnSvnLJF = ljnrgf
          OldSvnSvnCB = cnrgb
          OldSvnSvnLJB = ljnrgb
          OldSvnSvnDist = mindist
        else if (savenrg.eq.NEW) then
          NewSvnSvnC = cnrg
          NewSvnSvnLJ = ljnrg
          NewSvnSvnCF = cnrgf
          NewSvnSvnLJF = ljnrgf
          NewSvnSvnCB = cnrgb
          NewSvnSvnLJB = ljnrgb
          NewSvnSvnDist = mindist
        endif

c       now the protein-solvent energy
        if (NProteins.gt.0) then
          do i=1,NProteins
c           add the total solvent-protein energy for protein i
            nrg = proteinSolventEnergy(i,
     .              cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
            nrgf = nrgf + cnrgf + ljnrgf
            nrgb = nrgb + cnrgb + ljnrgb
            totalEnergy = totalEnergy + nrg
            if (savenrg.eq.OLD) then
              OldProSvnC(i) = cnrg
              OldProSvnLJ(i) = ljnrg
              OldProSvnCF(i) = cnrgf
              OldProSvnLJF(i) = ljnrgf
              OldProSvnCB(i) = cnrgb
              OldProSvnLJB(i) = ljnrgb
              OldProSvnDist(i) = mindist
            else if (savenrg.eq.NEW) then
              NewProSvnC(i) = cnrg
              NewProSvnLJ(i) = ljnrg
              NewProSvnCF(i) = cnrgf
              NewProSvnLJF(i) = ljnrgf
              NewProSvnCB(i) = cnrgb
              NewProSvnLJB(i) = ljnrgb
              NewProSvnDist(i) = mindist
            endif
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
              
              nrg = soluteSolventEnergy(DONTCARE,i,
     .              cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist,
     .              scl,sclf,sclb)
c             scale the energies
              nrg = nrg*scl
              cnrg = cnrg*scl
              ljnrg = ljnrg*scl
              cnrgf = cnrgf*sclf
              ljnrgf = ljnrgf*sclf
              cnrgb = cnrgb*sclb
              ljnrgb = ljnrgb*sclb
c             If the soft-core is on, the coulombic is scaled further here
              if (Solutesoftcore(i)) then
                cnrg = cnrg * (scl)**CoulSoftPower
                cnrgf = cnrgf * (sclf)**CoulSoftPower
                cnrgb = cnrgb * (sclb)**CoulSoftPower
              endif
              nrgf = nrgf + cnrgf + ljnrgf
              nrgb = nrgb + cnrgb + ljnrgb
c             totalEnergy = totalEnergy + nrg
              totalEnergy = totalEnergy + cnrg + ljnrg
              
              if (savenrg.eq.OLD) then
                OldSolSvnC(i) = cnrg
                OldSolSvnLJ(i) = ljnrg
                OldSolSvnCF(i) = cnrgf
                OldSolSvnLJF(i) = ljnrgf
                OldSolSvnCB(i) = cnrgb
                OldSolSvnLJB(i) = ljnrgb
                OldSolSvnDist(i) = mindist
                if (i.eq.SvnPrefSol .and. PrefSample) then
c                 save the solvent probabilities...
                  do j=1,NSolvents
                    SvnProb(j) = SvnNewProb(j)
                  enddo
                  SvnSumProb = SvnNewSumProb
                  SvnMaxProb = SvnNewMaxProb
                endif
              else if (savenrg.eq.NEW) then
                NewSolSvnC(i) = cnrg
                NewSolSvnLJ(i) = ljnrg
                NewSolSvnCF(i) = cnrgf
                NewSolSvnLJF(i) = ljnrgf
                NewSolSvnCB(i) = cnrgb
                NewSolSvnLJB(i) = ljnrgb
               NewSolSvnDist(i) = mindist
              endif      
            endif
          enddo
        endif
      endif

      if (savenrg.eq.OLD) then
        OldTotalEnergy = totalEnergy
        OldTotalEnergyF = nrgf
        OldTotalEnergyB = nrgb
      else if (savenrg.eq.NEW) then
        NewTotalEnergy = totalEnergy
        NewTotalEnergyF = nrgf
        NewTotalEnergyB = nrgb
      endif

c      write(6,*)totalEnergy,nrgf,nrgb

      if (.not.DifferentLambdas) then
c       if lambdas are all the same then reference and perturbed energies should be the same!      
        if (dabs(totalEnergy - nrgf).gt.dabs(TINY*totalEnergy) .or. 
     .                              dabs(totalEnergy - nrgb).gt.dabs(TINY*totalEnergy)) then
          write(printstring,*) "Different energies despite same lambdas! ",totalEnergy," vs. ",
     .                          nrgb," vs. ",nrgf
          call printLine(WARNING,printstring)
          call printLine(WARNING,"This is a serious bug in the program that means that the energy")
          call printLine(WARNING,"calculation is wrong. You should either try to debug the problem")
          call printLine(WARNING,"or try to submit a bug report.")
        endif
      endif
          
      return

      end
