      double precision function zeroGCsoluteEnergy(savenrg,igcs,cnrgf,cnrgb,ljnrgf,ljnrgb)
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
      
      ZeroGCSoluteEnergy= ZERO
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
c$$$      inrg = GCsoluteGCsoluteEnergy(savenrg,igcs,
c$$$     .               icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,mindist)
c$$$      ZeroGCSoluteEnergy= ZeroGCSoluteEnergy+ inrg
      inrg = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      
      
      if (savenrg.eq.NEW) then
        NewIndGCSGCSC = ZERO
        NewIndGCSGCSLJ = ZERO
        NewIndGCSGCSCF = ZERO
        NewIndGCSGCSLJF = ZERO
        NewIndGCSGCSCB = ZERO
        NewIndGCSGCSLJB = ZERO
      else if (savenrg.eq.OLD) then
        OldIndGCSGCSC = ZERO
        OldIndGCSGCSLJ = ZERO
        OldIndGCSGCSCF = ZERO
        OldIndGCSGCSLJF = ZERO
        OldIndGCSGCSCB = ZERO
        OldIndGCSGCSLJB = ZERO
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

c$$$          do j=1,nres
c$$$            ZeroGCSoluteEnergy= ZeroGCSoluteEnergy+ residueSingleGCsoluteEnergy(i,j,igcs,
c$$$     .            icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
c$$$            cnrg = cnrg + icnrg
c$$$            ljnrg = ljnrg + iljnrg
c$$$            cnrgf = cnrgf + icnrgf
c$$$            cnrgb = cnrgb + icnrgb
c$$$            ljnrgf = ljnrgf + iljnrgf
c$$$            ljnrgb = ljnrgb + iljnrgb
c$$$            iicnrg = iicnrg + icnrg
c$$$            iiljnrg = iiljnrg + iljnrg
c$$$            iicnrgf = iicnrgf + icnrgf
c$$$            iicnrgb = iicnrgb + icnrgb
c$$$            iiljnrgf = iiljnrgf + iljnrgf
c$$$            iiljnrgb = iiljnrgb + iljnrgb
c$$$            if (imindist.lt.mindist) mindist = imindist
c$$$          enddo
c           PRINT*, "ZeroGCSoluteEnergyafter res" ,GCsoluteEnergy
          if (savenrg.eq.NEW) then
            NewIndGCSProC(i) = ZERO
            NewIndGCSProLJ(i) = ZERO
            NewIndGCSProCF(i) = ZERO
            NewIndGCSProLJF(i) = ZERO
            NewIndGCSProCB(i) = ZERO
            NewIndGCSProLJB(i) = ZERO
          else if (savenrg.eq.OLD) then
            OldIndGCSProC(i) = ZERO
c            PRINT*, "OldIndGCSProC in gcsolnrg" ,OldIndGCSProC(i)
            OldIndGCSProLJ(i) = ZERO
            OldIndGCSProCF(i) = ZERO
            OldIndGCSProLJF(i) = ZERO
            OldIndGCSProCB(i) = ZERO
            OldIndGCSProLJB(i) = ZERO
          endif
        enddo
      endif
         
      
c       This looks at solvent-solute energies
       if (NSolvents.gt.0) then
        do i=1,NSolvents
c$$$        ZeroGCSoluteEnergy= ZeroGCSoluteEnergy+ GCsoluteSingleSolventEnergy(igcs,i,
c$$$     .                    icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
c$$$             
c$$$            cnrg = cnrg + icnrg
c$$$            ljnrg = ljnrg + iljnrg
c$$$            cnrgf = cnrgf + icnrgf
c$$$            cnrgb = cnrgb + icnrgb
c$$$            ljnrgf = ljnrgf + iljnrgf
c$$$            ljnrgb = ljnrgb + iljnrgb

        if (savenrg.eq.OLD) then
          OldIndGCSSvnC(i) = ZERO
          OldIndGCSSvnLJ(i) = ZERO
          OldIndGCSSvnCF(i) = ZERO
          OldIndGCSSvnLJF(i) = ZERO
          OldIndGCSSvnCB(i) = ZERO
          OldIndGCSSvnLJB(i) = ZERO
          OldGCSSvnDist(i) = ZERO
        else if (savenrg.eq.NEW) then
          NewIndGCSSvnC(i) = ZERO
          NewIndGCSSvnLJ(i) = ZERO
          NewIndGCSSvnCF(i) = ZERO
          NewIndGCSSvnLJF(i) = ZERO
          NewIndGCSSvnCB(i) = ZERO
          NewIndGCSSvnLJB(i) = ZERO
          NewGCSSvnDist(i) = ZERO
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

c$$$          nrg = soluteSingleGCSoluteEnergy(i,igcs,
c$$$     .        icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist,
c$$$     .        scl,sclf,sclb)
          nrg = ZERO
          if (Solutesoftcore(i)) then
            icnrg = icnrg * (scl)**CoulSoftPower
            icnrgf = icnrgf * (sclf)**CoulSoftPower
            icnrgb = icnrgb * (sclb)**CoulSoftPower
          endif
  
c$$$          cnrg = cnrg + icnrg*scl
c$$$          ljnrg = ljnrg + iljnrg*scl
c$$$          cnrgf = cnrgf + icnrgf*sclf
c$$$          cnrgb = cnrgb + icnrgb*sclb
c$$$          ljnrgf = ljnrgf + iljnrgf*sclf
c$$$          ljnrgb = ljnrgb + iljnrgb*sclb
c$$$          ZeroGCSoluteEnergy= ZeroGCSoluteEnergy+ scl*icnrg + scl*iljnrg
          if (savenrg.eq.OLD) then
            OldIndGCSSolC(i) = ZERO
            OldIndGCSSolLJ(i) = ZERO
            OldIndGCSSolCF(i) = ZERO
            OldIndGCSSolLJF(i) = ZERO
            OldIndGCSSolCB(i) = ZERO
            OldIndGCSSolLJB(i) = ZERO
          else if (savenrg.eq.NEW) then
            NewIndGCSSolC(i) = ZERO
            NewIndGCSSolLJ(i) = ZERO
            NewIndGCSSolCF(i) = ZERO
            NewIndGCSSolLJF(i) = ZERO
            NewIndGCSSolCB(i) = ZERO
            NewIndGCSSolLJB(i) = ZERO
          endif
        enddo
      endif

      if (savenrg.eq.NEW) then
        NewIndGCSC = ZERO
c        PRINT*, "NewIndGCSC" ,NewIndGCSC
        NewIndGCSLJ = ZERO
        NewIndGCSCF = ZERO
        NewIndGCSLJF = ZERO
        NewIndGCSCB = ZERO
        NewIndGCSLJB = ZERO
      else if (savenrg.eq.OLD) then
        OldIndGCSC = ZERO
c        PRINT*, "OldIndGCSC" ,OldIndGCSC
        OldIndGCSLJ = ZERO
        OldIndGCSCF = ZERO
        OldIndGCSLJF = ZERO
        OldIndGCSCB = ZERO
        OldIndGCSLJB = ZERO
      endif
      nrgf = ZERO
      nrgb = ZERO
      
      return
      end
