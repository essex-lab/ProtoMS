      double precision function TotalSurface(savenrg,surfacef,surfaceb,
     .  sanrg,sanrgf,sanrgb)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'surface.inc'
      include 'printer.inc'
      include 'simparams.inc'
      include 'enums.inc'
      include 'flags.inc'
c################################################
c     Calculate the Surface Area of the system 
c     
c     Julien Michel November 2004
c################################################      
      integer savenrg
      double precision temp,surfacef,surfaceb
      double precision sanrg,sanrgf,sanrgb
      integer i,j,nres
      double precision SoluteSurface
      double precision ResidueSurface
      double precision SolventSurface
      double precision solsurfacef,solsurfaceb,solsanrg,solsanrgf,solsanrgb
      double precision svnsurfacef,svnsurfaceb,svnsanrg,svnsanrgf,svnsanrgb
      double precision ressurfacef,ressurfaceb,ressanrg,ressanrgf,ressanrgb
c      double precision scl,sclf,sclb
c      double precision getLowScale,getHighScale

      TotalSurface = OldTotalSurface
      surfacef = OldTotalSurfaceF
      surfaceb = OldTotalSurfaceB
      sanrg = oldSAEnergy
      sanrgf = oldSAEnergyF
      sanrgb = oldSAEnergyB

c     If an OLD surface is requested, then it is already 
c     known and can be immediately returned
      if (savenrg .eq. OLD) return

c     Mark which residues needs to have their surface area updated
c     

      if (.not. DualredoSA) call markSAupdate
      
      temp = ZERO
      if (NProteins .gt. ZERO) then      
        do 100 i=1,NProteins
          nres = NResidues(i)
          do 101 j=1,nres
            if (ResidueredoSA(i,j)) then
              temp = ResidueSurface(savenrg,i,j,ressurfacef,ressurfaceb,
     .          ressanrg, ressanrgf,ressanrgb)
              TotalSurface = TotalSurface + temp
              surfacef = surfacef + ressurfacef
              surfaceb = surfaceb + ressurfaceb
c             write(printstring,*) ' ressa ',ressanrg,' f ',ressanrgf,
c     .          ' b' ,ressanrgb,' surf ',temp,' f ',solsurfacef,' b ',solsurfaceb
c             call printLine(DEBUG,printstring)
              sanrg = sanrg + ressanrg
              sanrgf = sanrgf + ressanrgf
              sanrgb = sanrgb + ressanrgb
              if (savenrg .eq. OLD) then
                OldResidueSurf(i,j) = temp
                OldResidueSurfF(i,j) = ressurfacef
                OldResidueSurfB(i,j) = ressurfaceb
              else if (savenrg .eq. NEW) then
                NewResidueSurf(i,j) = temp
                NewResidueSurfF(i,j) = ressurfacef
                NewResidueSurfB(i,j) = ressurfaceb                
              endif
          endif
 101      continue
 100    continue
      endif
c      write(printstring,*) 'proteins  surf ',Totalsurface
c      call printLine(DEBUG,printstring)
      if (NSolutes .gt. ZERO) then
        do 200 i=1,NSolutes
          if (IgnoreSolute(i)) goto 200
          if (SoluteredoSA(i)) then
            temp = SoluteSurface(savenrg,i,solsurfacef,solsurfaceb,
     .        solsanrg,solsanrgf,solsanrgb) 
c            write(printstring,*)  ' solute ',i,' surf ',temp,' f ' ,solsurfacef,' b ',solsurfaceb
c            call printLine(DEBUG,printstring)
            TotalSurface = TotalSurface + temp 
            surfacef = surfacef + solsurfacef 
            surfaceb = surfaceb + solsurfaceb 
            sanrg = sanrg + solsanrg  
            sanrgf = sanrgf + solsanrgf 
            sanrgb = sanrgb + solsanrgb 
c           write(printstring,*) ' solsanrg ',solsanrg,' f ',solsanrgf,
c     .        ' b' ,solsanrgb,' surf ',temp,' f ',solsurfacef,' b ',solsurfaceb
c           call printLine(DEBUG,printstring)
            if (savenrg .eq. OLD) then
              OldSoluteSurf(i) = temp 
              OldSoluteSurfF(i) = solsurfacef
              OldSoluteSurfB(i) = solsurfaceb
            else if (savenrg .eq. NEW) then
              NewSoluteSurf(i) = temp
              NewSoluteSurfF(i) = solsurfacef
              NewSoluteSurfB(i) = solsurfaceb                
            endif
          endif
 200    continue
      endif

      if (NSolvents .gt. ZERO) then
        do 300 i=1,NSolvents
          if (SolventredoSA(i)) then
            temp = SolventSurface(savenrg,i,svnsurfacef,svnsurfaceb,
     .        svnsanrg,svnsanrgf,svnsanrgb)
c            write(printstring,*) 'svn ',i,' surf ',temp
c            call printLine(DEBUG,printstring)
            TotalSurface = TotalSurface + temp
            surfacef = surfacef + svnsurfacef
            surfaceb = surfaceb + svnsurfaceb
            sanrg = sanrg + svnsanrg
            sanrgf = sanrgf + svnsanrgf 
            sanrgb = sanrgb + svnsanrgb
c           write(printstring,*) ' svnsanrg ',svnsanrg,' f ',svnsanrgf,
c    .        ' b' ,svnsanrgb,' surf ',temp,' f ',svnsurfacef,' b ',svnsurfaceb
c           call printLine(DEBUG,printstring)
            if (savenrg .eq. OLD) then
              OldSolventSurf(i) = temp
              OldSolventSurfF(i) = svnsurfacef
              OldSolventSurfB(i) = svnsurfaceb
            else if (savenrg .eq. NEW) then
              NewSolventSurf(i) = temp
              NewSolventSurfF(i) = svnsurfacef
              NewSolventSurfB(i) = svnsurfaceb                
              endif
          endif
 300    continue
      endif

c     Need to save the surface as well as the surface area energy
      if (savenrg .eq. OLD) then
        oldTotalSurface = TotalSurface
        oldTotalSurfaceF = surfacef
        oldTotalSurfaceB = surfaceb
        oldSAEnergy = sanrg
        oldSAEnergyF = sanrgf
        oldSAEnergyB = sanrgb
      else if (savenrg .eq. NEW) then
        NewTotalSurface = TotalSurface
        NewTotalSurfaceF = surfacef
        NewTotalSurfaceB = surfaceb
        NewSAEnergy = sanrg
        NewSAEnergyF = sanrgf
        NewSAEnergyB = sanrgb
      endif
      

      return
      end
