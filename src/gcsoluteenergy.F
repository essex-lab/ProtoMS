      double precision function GCsoluteEnergy(savenrg,igcs,cnrgf,cnrgb,ljnrgf,ljnrgb)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'flags.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'solvents.inc'
      include 'simparams.inc'
      include 'enums.inc'
      include 'energies.inc'
c##############################################################
c
c     This function returns the total energy of GC-Solute igcs,
c     (bearing in mind what has come before)
c
c     If savenrg = 1 then this is the old energy, else if it 
c     is 2, then it is the new energy, else don't save the
c     energy!
c
c     (C) MSB 22/04/09
c
c#################################################################

      integer i,j,savenrg,igcs,nres
      double precision cnrg,ljnrg,mindist,icnrg,iljnrg,imindist,inrg
      double precision iicnrg,iiljnrg,iicnrgf,iicnrgb,iiljnrgf,iiljnrgb
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb,icnrgf,icnrgb,iljnrgf,iljnrgb
      double precision nrg,scl,sclf,sclb
      double precision soluteSingleGCSoluteEnergy
      double precision GCsoluteGCsoluteEnergy,GCsoluteSingleSolventEnergy
      double precision residueSingleGCsoluteEnergy,getHighScale,getLowScale
      double precision nrgb,nrgf
      
      GCsoluteEnergy = ZERO
      nrgf = ZERO
      nrgb = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      ljnrgf = ZERO
      cnrgb = ZERO
      ljnrgb = ZERO
      
      if (NGCSolutes.le.0 .or. igcs.gt.MAXSOLVENTS) return
      
      if (DoneGCSNrg(igcs)) then
        write(printstring,*) "Have we already done this solute? - ",igcs
        call printLine(WARNING,printstring)
        return
      endif
      
c     This is the bit which works out solute-solute energies
      inrg = GCsoluteGCsoluteEnergy(savenrg,igcs,
     .               icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,mindist)
      GCsoluteEnergy = GCsoluteEnergy + inrg
      cnrg = icnrg
      ljnrg = iljnrg
      cnrgf = icnrgf
      cnrgb = icnrgb
      ljnrgf = iljnrgf
      ljnrgb = iljnrgb
      
      
      if (savenrg.eq.NEW) then
        NewIndGCSGCSC = icnrg
        NewIndGCSGCSLJ = iljnrg
        NewIndGCSGCSCF = icnrgf
        NewIndGCSGCSLJF = iljnrgf
        NewIndGCSGCSCB = icnrgb
        NewIndGCSGCSLJB = iljnrgb
      else if (savenrg.eq.OLD) then
        OldIndGCSGCSC = icnrg
        OldIndGCSGCSLJ = iljnrg
        OldIndGCSGCSCF = icnrgf
        OldIndGCSGCSLJF = iljnrgf
        OldIndGCSGCSCB = icnrgb
        OldIndGCSGCSLJB = iljnrgb
      endif
      
c      Debugging output
      
c     Now, it looks at protein-solute energies
      if (NProteins.gt.0) then
         do i=1,NProteins
          nres = NResidues(i)
          iicnrg = ZERO
          iiljnrg = ZERO
          iicnrgf = ZERO
          iiljnrgf = ZERO
          iicnrgb = ZERO
          iiljnrgb = ZERO

          do j=1,nres
            GCsoluteEnergy = GCsoluteEnergy + residueSingleGCsoluteEnergy(i,j,igcs,
     .            icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
            cnrg = cnrg + icnrg
            ljnrg = ljnrg + iljnrg
            cnrgf = cnrgf + icnrgf
            cnrgb = cnrgb + icnrgb
            ljnrgf = ljnrgf + iljnrgf
            ljnrgb = ljnrgb + iljnrgb
            iicnrg = iicnrg + icnrg
            iiljnrg = iiljnrg + iljnrg
            iicnrgf = iicnrgf + icnrgf
            iicnrgb = iicnrgb + icnrgb
            iiljnrgf = iiljnrgf + iljnrgf
            iiljnrgb = iiljnrgb + iljnrgb
            if (imindist.lt.mindist) mindist = imindist
          enddo
c           PRINT*, "GCsoluteEnergy after res" ,GCsoluteEnergy
          if (savenrg.eq.NEW) then
            NewIndGCSProC(i) = iicnrg         
            NewIndGCSProLJ(i) = iiljnrg
            NewIndGCSProCF(i) = iicnrgf
            NewIndGCSProLJF(i) = iiljnrgf
            NewIndGCSProCB(i) = iicnrgb
            NewIndGCSProLJB(i) = iiljnrgb
          else if (savenrg.eq.OLD) then
            OldIndGCSProC(i) = iicnrg
c            PRINT*, "OldIndGCSProC in gcsolnrg" ,OldIndGCSProC(i)
            OldIndGCSProLJ(i) = iiljnrg
            OldIndGCSProCF(i) = iicnrgf
            OldIndGCSProLJF(i) = iiljnrgf
            OldIndGCSProCB(i) = iicnrgb
            OldIndGCSProLJB(i) = iiljnrgb
          endif
        enddo
      endif
         
      
c       This looks at solvent-solute energies

       if (NSolvents.gt.0) then
        do i=1,NSolvents
        GCSoluteEnergy = GCSoluteEnergy + GCsoluteSingleSolventEnergy(igcs,i,
     .                    icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
             
            cnrg = cnrg + icnrg
            ljnrg = ljnrg + iljnrg
            cnrgf = cnrgf + icnrgf
            cnrgb = cnrgb + icnrgb
            ljnrgf = ljnrgf + iljnrgf
            ljnrgb = ljnrgb + iljnrgb

        if (savenrg.eq.OLD) then
          OldIndGCSSvnC(i) = icnrg
          OldIndGCSSvnLJ(i) = iljnrg
          OldIndGCSSvnCF(i) = icnrgf
          OldIndGCSSvnLJF(i) = iljnrgf
          OldIndGCSSvnCB(i) = icnrgb
          OldIndGCSSvnLJB(i) = iljnrgb
          OldGCSSvnDist(i) = imindist
        else if (savenrg.eq.NEW) then
          NewIndGCSSvnC(i) = icnrg
          NewIndGCSSvnLJ(i) = iljnrg
          NewIndGCSSvnCF(i) = icnrgf
          NewIndGCSSvnLJF(i) = iljnrgf
          NewIndGCSSvnCB(i) = icnrgb
          NewIndGCSSvnLJB(i) = iljnrgb
          NewGCSSvnDist(i) = imindist
        endif
       enddo
       endif     
c       PRINT*, "nrg in gcsolute" ,GCsoluteEnergy,cnrgf+ljnrgf,cnrgb+ljnrgb

      if (NSolutes.gt.0) then
        do i=1,NSolutes

c         see whether we need to scale the energies to account for 
c         this being a dual topology solute
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

          nrg = soluteSingleGCSoluteEnergy(i,igcs,
     .        icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist,
     .        scl,sclf,sclb)

          if (Solutesoftcore(i)) then
            icnrg = icnrg * (scl)**CoulSoftPower
            icnrgf = icnrgf * (sclf)**CoulSoftPower
            icnrgb = icnrgb * (sclb)**CoulSoftPower
          endif
  
          cnrg = cnrg + icnrg*scl
          ljnrg = ljnrg + iljnrg*scl
          cnrgf = cnrgf + icnrgf*sclf
          cnrgb = cnrgb + icnrgb*sclb
          ljnrgf = ljnrgf + iljnrgf*sclf
          ljnrgb = ljnrgb + iljnrgb*sclb
          GCSoluteEnergy = GCSoluteEnergy + scl*icnrg + scl*iljnrg
          if (savenrg.eq.OLD) then
            OldIndGCSSolC(i) = icnrg*scl
            OldIndGCSSolLJ(i) = iljnrg*scl
            OldIndGCSSolCF(i) = icnrgf*sclf
            OldIndGCSSolLJF(i) = iljnrgf*sclf
            OldIndGCSSolCB(i) = icnrgb*sclb
            OldIndGCSSolLJB(i) = iljnrgb*sclb
          else if (savenrg.eq.NEW) then
            NewIndGCSSolC(i) = icnrg*scl
            NewIndGCSSolLJ(i) = iljnrg*scl
            NewIndGCSSolCF(i) = icnrgf*sclf
            NewIndGCSSolLJF(i) = iljnrgf*sclf
            NewIndGCSSolCB(i) = icnrgb*sclb
            NewIndGCSSolLJB(i) = iljnrgb*sclb
            NewGCSSvnDist(i) = imindist
          endif
        enddo
      endif

      if (savenrg.eq.NEW) then
        NewIndGCSC = cnrg
c        PRINT*, "NewIndGCSC" ,NewIndGCSC
        NewIndGCSLJ = ljnrg        
        NewIndGCSCF = cnrgf
        NewIndGCSLJF = ljnrgf
        NewIndGCSCB = cnrgb
        NewIndGCSLJB = ljnrgb
      else if (savenrg.eq.OLD) then
        OldIndGCSC = cnrg
c        PRINT*, "OldIndGCSC" ,OldIndGCSC
        OldIndGCSLJ = ljnrg  
        OldIndGCSCF = cnrgf
        OldIndGCSLJF = ljnrgf
        OldIndGCSCB = cnrgb
        OldIndGCSLJB = ljnrgb
      endif
      nrgf = cnrgf + ljnrgf
      nrgb = cnrgb + ljnrgb
      
      return
      end
