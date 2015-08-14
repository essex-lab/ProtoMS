      subroutine dualtopologysasa(isol,savenrg,sanrg,sanrgf,sanrgb,surface,surfacef,surfaceb)
      include 'dimensions.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'constants.inc'
      include 'printer.inc'
      include 'extraenergies.inc'
      include 'solutes.inc'
      include 'simparams.inc'
c##################################################
c     Allow correct SASA calculation for dual 
c     topology simulations. 
c     Note that the following only works 
c     for ONE pair of dual topology solutes
c
c
c##################################################      
      integer isol,jsol
      integer savenrg
      double precision TotalSurface
      double precision temp,surfacef,surfaceb
      double precision sanrg,sanrgf,sanrgb
      double precision temp1,surfacef1,surfaceb1
      double precision sanrg1,sanrgf1,sanrgb1
      double precision temp2,surfacef2,surfaceb2
      double precision sanrg2,sanrgf2,sanrgb2
      double precision getLowScale,getHighScale
      double precision surface
      double precision lowscale,lowscalef,lowscaleb
      double precision highscale,highscalef,highscaleb      

      double precision svOldTotalsurface,svOldTotalSurfaceF,svOldTotalSurfaceB
      double precision svsanrg,svsanrgf,svsanrgb

      surface = OldTotalSurface
      surfacef = OldTotalSurfaceF
      surfaceb = OldTotalSurfaceB
      sanrg = oldSAEnergy
      sanrgf = oldSAEnergyF
      sanrgb = oldSAEnergyB 

      if (savenrg .eq. OLD) return

      DualredoSA = .true.

      jsol = DualTopologySolute(isol)
c
c     Save Old values
c
      svOldTotalSurface = OldTotalSurface
      svOldTotalSurfaceF = OldTotalSurfaceF
      svOldTotalSurfaceB = OldTotalSurfaceB
      svsanrg = oldSAEnergy
      svsanrgf = oldSAEnergyF
      svsanrgb = oldSAEnergyB

c     turn off solutes
      IgnoreSolute(jsol) = .true.
      IgnoreSolute(isol) = .false.
c
c     initialise the SA energy
c
      call redosasa

      temp = TotalSurface(savenrg,surfacef,surfaceb,sanrg,sanrgf,sanrgb)
c      write(printstring,*) ' SASA ignore ',jsol,' surf ',temp
c      call printLine(DEBUG,printstring) 

c     Save sanrg/sanrgf/sanrgb temp surfacef surfaceb
      temp1     = temp
      surfacef1 = surfacef
      surfaceb1 = surfaceb
      sanrg1    = sanrg
      sanrgf1   = sanrgf
      sanrgb1   = sanrgb

c     turn back  DualTopologySolute(i) on and i off
      IgnoreSolute(jsol) = .false.            
      IgnoreSolute(isol) = .true.
c
c     initialise the SA energy
c
      call redosasa

      temp = TotalSurface(savenrg,surfacef,surfaceb,sanrg,sanrgf,sanrgb)
c      write(printstring,*) ' SASA ignore ',isol,' surf ',temp
c      call printLine(DEBUG,printstring) 
c     Save sanrg/sanrgf/sanrgb temp surfacef surfaceb
      temp2     = temp
      surfacef2 = surfacef
      surfaceb2 = surfaceb
      sanrg2    = sanrg
      sanrgf2   = sanrgf
      sanrgb2   = sanrgb
c     turn back i on 
      IgnoreSolute(jsol) = .false.            
      IgnoreSolute(isol) = .false.
c     ***********************
c     Calculate actual SASA values from linear combination 
      lowscale = getLowScale(Lambda) 
      lowscalef = getLowScale(LambdaF) 
      lowscaleb = getLowScale(LambdaB) 
      highscale = getHighScale(Lambda)
      highscalef = getHighScale(LambdaF)
      highscaleb = getHighScale(LambdaB)

      surface = temp1 * lowscale + temp2 * highscale
      surfacef =  surfacef1 * lowscalef +  surfacef2 * highscalef
      surfaceb =  surfaceb1 * lowscaleb + surfaceb2 * highscaleb
      sanrg =  sanrg1 * lowscale + sanrg2 * highscale
      sanrgf = sanrgf1 * lowscalef + sanrgf2 * highscalef
      sanrgb = sanrgb1 * lowscaleb + sanrgb2 * highscaleb
c
c     Regenerate Old values for complete system
c
      OldTotalSurface = svOldTotalSurface 
      OldTotalSurfaceF = svOldTotalSurfaceF
      OldTotalSurfaceB = svOldTotalSurfaceB
      oldSAEnergy = svsanrg
      oldSAEnergyF = svsanrgf
      oldSAEnergyB = svsanrgb
c     THEN SAVE IN THE NEW SECTION 
c     Need to save the surface as well as the surface area energy
      NewTotalSurface = surface
      NewTotalSurfaceF = surfacef
      NewTotalSurfaceB = surfaceb
      NewSAEnergy = sanrg
      NewSAEnergyF = sanrgf
      NewSAEnergyB = sanrgb

c      write(printstring,*) ' SolSurfidual ',NewDualSoluteSurf(idual),
c     .  ' atom 1 ',NewDualSolAtomSurf(idual,1),
c     .  ' SolSurfJdual ',NewDualSoluteSurf(jdual),
c     .  ' atom 1 ',NewDualSolAtomSurf(jdual,1),
c     .  ' SolSurf ',NewSoluteSurf(isol)
c      call printLine(DEBUG,printstring)

c      STOP
      
      DualredoSA = .false.

      return
      end


      subroutine redosasa
      include 'dimensions.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'solvents.inc'
      include 'energies.inc'
      include 'flags.inc'
      include 'constants.inc'
      include 'printer.inc'
      include 'surface.inc'     

      integer i,j,k
      integer nats,nres,resstart,resend

      oldTotalSurface = ZERO
      oldTotalSurfaceF = ZERO
      oldTotalSurfaceB = ZERO
      oldSAEnergy = ZERO
      oldSAEnergyF = ZERO
      oldSAEnergyB = ZERO
      NewTotalSurface = ZERO
      NewTotalSurfaceF = ZERO
      NewTotalSurfaceB = ZERO
      NewSAEnergy = ZERO
      NewSAEnergyF = ZERO
      NewSAEnergyB = ZERO

      if (doSA) then 
        call clearSurfaces
      endif

      if (NSolutes .gt. 0) then
        do 1 i=1,NSolutes
          if (IgnoreSolute(i)) goto 1
          SoluteredoSA(i) = .true.
          nats = NSoLAtoms(i)
          do 2 j=i,nats
            SolAtredoSA(i,j) = .true.
 2        continue
 1     continue
      endif

      if (NSolvents .gt. 0) then
        do 3 i=1,NSolvents
          SolventredoSA(i) = .true.
          nats = NSvnAtoms(i)
          do 4 j=1,nats
            SvnAtredoSA(i,j) = .true.            
 4        continue
 3      continue
      endif

      if (NProteins .gt. 0) then
        do 5 i=1,NProteins
          nres = NResidues(i)
          do 6 j=1,nres
            ResidueredoSA(i,j) = .true.
            resstart = BBZid(i,j,1)
            resend = BBZid(i,j,4)            
            do 7 k=resstart,resend
              ProtAtredoSA(i,k) = .true.
 7          continue
 6        continue
 5      continue
      endif     
      
      return
      end
