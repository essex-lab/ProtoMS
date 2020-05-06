      double precision function SolventSurface(savenrg,isvn,surfacef,surfaceb,sanrg,sanrgf,sanrgb)
c
c     Compute the Surface Area of one solute
c
      include 'dimensions.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'printer.inc'
      include 'surface.inc'
      include 'flags.inc'
      include 'parameters.inc'
      include 'enums.inc'
      include 'simparams.inc'
      integer savenrg
      integer isvn,ipar
      integer i,j,k,l,m,nats,natsother,nres,resstart,resend
      double precision surfacef,surfaceb,diffsurfF,diffsurfB
      double precision sanrg,sanrgf,sanrgb
      double precision ipoint(nsurfpoints*3)
      double precision iarea(nsurfpoints),darea
      logical  iskip(nsurfpoints)
      double precision ipointf(nsurfpoints*3),iareaf(nsurfpoints),dareaf
      logical iskipf(nsurfpoints)
      double precision ipointb(nsurfpoints*3),iareab(nsurfpoints),dareab
      logical iskipb(nsurfpoints)
      double precision ri,ri2,rj,rj2,dij2,dkj2
      double precision xi,yi,zi,xj,yj,zj
      double precision isurf,oldsurf,diffsurf
      double precision isurftens,isurftensf,isurftensb
      double precision isurff,isurfb
      double precision oldsurff,oldsurfb
      double precision rjf,rjf2,rjb,rjb2,dijf2,dijb2
      double precision rif,rif2,rib,rib2
      double precision xif,yif,zif,xib,yib,zib
      double precision xjf,yjf,zjf,xjb,yjb,zjb
      integer jpar
      logical jperturbedatom,jvanish,jvanishF,jvanishB

      SolventSurface = ZERO      
      surfacef = ZERO
      surfaceb = ZERO
      sanrg = ZERO
      sanrgf = ZERO
      sanrgb = ZERO

      nats = NSvnAtoms(isvn)

      do 10 i=1,nats
c       For each atom i
c       Check if it is necessary to compute the surface again
        if (.not. SvnAtredoSA(isvn,i)) then
c          isurf = OldSvnAtomSurf(isvn,i)
c          write(printstring,*) ' svn ',isvn,' atom ',i,
c     &      'SA did not change, regenerated ',isurf
c          call printLine(DEBUG,printstring)
c          goto 201          
          goto 10
        endif
c       Skip atom if radii is 0
c        ri = SvnAtRadius(isvn,i)
        ipar = ZSvnPar(isvn,i)
        if (ipar .lt. 0) then
c         Perturbed atom
          ri = VariableSurfRadius(-ipar,2)
          isurftens = VariableSurfTension(-ipar,2)
          write(printstring,*) 'No support for perturbed Solvent atom'
          call printLine(FATAL,printstring)
          call closeProgram(-1)
        endif
        ri = SurfRadDB(CLJType(ipar))
        if (ri .lt. 0.01d0) goto 10
        isurftens = SurfTensionDB(CLJType(ipar))
        isurftensf = isurftens
        isurftensb = isurftens
        ri = ri + probe
        ri2 = ri **2
        xi = ZSvnCoords(isvn,i,1)
        yi = ZSvnCoords(isvn,i,2)
        zi = ZSvnCoords(isvn,i,3)
c        write(printstring,*) 'isurf ',isurftens,' f ',isurftensf,' b ',isurftensb
c        call printLine(DEBUG,printstring)
        darea = ( FOURPI * ri2 ) / nsurfpoints
        xif = xi
        yif = yi
        zif = zi
        rif = ri
        rif2 = ri2
        dareaf = darea
        xib = xi
        yib = yi
        zib = zi
        rib = ri
        rib2 = ri2
        dareab = darea
c       Set up the surface elements of atom i
        if (SAquality .eq. 1) then
          do 20 j=1,nsurfpoints
            ipoint(3*j-2) = points12(3*j-2) * ri + xi
            ipoint(3*j-1) = points12(3*j-1) * ri + yi
            ipoint(3*j) = points12(3*j) * ri + zi
            iarea(j) = darea
            iskip(j) = .false.
 20       continue
        else if (SAquality .eq. 2) then
          do 21 j=1,nsurfpoints
            ipoint(3*j-2) = points42(3*j-2) * ri + xi
            ipoint(3*j-1) = points42(3*j-1) * ri + yi
            ipoint(3*j) = points42(3*j) * ri + zi
            iarea(j) = darea
            iskip(j) = .false.
 21       continue          
        else if (SAquality .eq. 3) then
          do 22 j=1,nsurfpoints
            ipoint(3*j-2) = points162(3*j-2) * ri + xi
            ipoint(3*j-1) = points162(3*j-1) * ri + yi
            ipoint(3*j) = points162(3*j) * ri + zi
            iarea(j) = darea
            iskip(j) = .false.
 22       continue          
        else if (SAquality .eq. 4) then
          do 23 j=1,nsurfpoints
            ipoint(3*j-2) = points642(3*j-2) * ri + xi
            ipoint(3*j-1) = points642(3*j-1) * ri + yi
            ipoint(3*j) = points642(3*j) * ri + zi
            iarea(j) = darea
            iskip(j) = .false.
 23       continue          
        endif
c       Plus forward/backward surface elements
        do 201 j=1,nsurfpoints
c         forward
          ipointf(3*j-2) = ipoint(3*j-2) 
          ipointf(3*j-1) = ipoint(3*j-1) 
          ipointf(3*j) = ipoint(3*j) 
          iareaf(j) = dareaf
          iskipf(j) = .false.
c         backward
          ipointb(3*j-2) = ipoint(3*j-2) 
          ipointb(3*j-1) = ipoint(3*j-1) 
          ipointb(3*j) = ipoint(3*j) 
          iareab(j) = dareab
          iskipb(j) = .false.
 201    continue
c       Now do Solvent with other Proteins
        if (NProteins .gt. ZERO) then
          do 60 j=1,NProteins
            nres = NResidues(j)
            do 70 k=1,nres
              resstart = BBZid(j,k,1)
              resend = BBZid(j,k,4)
              do 80 l=resstart,resend
                xj = ZCoords(j,l,1)
                yj = ZCoords(j,l,2)
                zj = ZCoords(j,l,3)
c                rj = ProtAtRadius(j,l)
                rj = SurfRadDB(CLJType(ZPar(j,l)))
c               Skip if radius is zero
                if (rj .lt. 0.01d0) goto 80
                rj = rj + probe
                rj2 = rj**2
                dij2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c               Skip atom if no overlap
                if (dij2 .gt. (ri+rj)**2) goto 80
c               Eliminate buried elements
                do 90 m=1,nsurfpoints
c                 skip element if already buried (by another atom then)
                  if (iskip(m)) goto 90
                  dkj2 = ( ipoint(3*m-2)-xj )**2 + ( ipoint(3*m-1)-yj )**2 +
     &              ( ipoint(3*m) - zj)**2
                  if (dkj2 .lt. rj2) then 
                    iskip(m) = .true.
                    iskipf(m) = .true.
                    iskipb(m) = .true.
                  endif
 90             continue
 80           continue
 70         continue
 60       continue
        endif
c       Now do Solvent with other solvents
          do 100 j=1,NSolvents
            natsother = NSvnAtoms(j)
            do 110 k=1,natsother
              if (i .eq. k .and. isvn .eq. j) goto 110
              xj = ZSvnCoords(j,k,1)
              yj = ZSvnCoords(j,k,2)
              zj = ZSvnCoords(j,k,3)
c              rj = SvnAtRadius(j,k)
              rj = SurfRadDB(CLJType(ZSvnPar(j,k)))
c             Skip if radius is zero
              if (rj .lt. 0.01d0) goto 110
              rj = rj + probe
              rj2 = rj**2
              dij2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c             Skip atom if no overlap
              if (dij2 .gt. (ri+rj)**2) goto 110
c             Eliminate buried elements
              do 120 l=1,nsurfpoints
c               skip element if already buried (by another atom then)
                if (iskip(l)) goto 120
                dkj2 = ( ipoint(3*l-2)-xj )**2 + ( ipoint(3*l-1)-yj )**2 +
     &            ( ipoint(3*l) - zj)**2
                if (dkj2 .lt. rj2) then 
                  iskip(l) = .true.
                  iskipf(l) = .true.
                  iskipb(l) = .true.
                endif
 120          continue               
 110        continue
 100      continue
c       Now do Solvent with other Solutes
        if (NSolutes .gt. ZERO) then
          do 30 j=1,NSolutes
            natsother = NSolAtoms(j)
            do 40 k=1,natsother
              jperturbedatom = .false.
              jvanish = .false.
              jvanishF = .false.
              jvanishB = .false.
              jpar = ZSolPar(j,k)
              if (jpar .lt. 0) then
                jperturbedatom = .true.
c               Perturbed Atom
                rj = VariableSurfRadius(-jpar,2)
                if (DifferentLambdas) then
c                 Set up forward/backward ri and surface tensions
                  rjf = VariableSurfRadius(-jpar,5)
                  rjb = VariableSurfRadius(-jpar,4)
                else
                  rjf = rj
                  rjb = rj 
                endif
              else
c               Normal Atom
                rj = SurfRadDB(CLJType(jpar))
                rjf = rj
                rjb = rj
              endif

              if (rj .lt. 0.01d0  ) jvanish = .true.
              if (rjf .lt. 0.01d0 ) jvanishF = .true.
              if (rjb .lt. 0.01d0 ) jvanishB = .true.
              rj = rj + probe
              rj2 = rj **2
              rjf = rjf + probe
              rjf2 = rjf ** 2
              rjb = rjb + probe
              rjb2 = rjb ** 2
              if (DualTopologySolute(j) .ne. 0) jperturbedatom = .true.
              xj = ZSolCoords(j,k,1)
              yj = ZSolCoords(j,k,2)
              zj = ZSolCoords(j,k,3)
              if (NSolVariable(j) .gt. 0) then
                jperturbedatom = .true.
                xjf = ZSolCoords(NSolForwards(j),k,1)
                yjf = ZSolCoords(NSolForwards(j),k,2)
                zjf = ZSolCoords(NSolForwards(j),k,3)
                xjb = ZSolCoords(NSolBackwards(j),k,1)
                yjb = ZSolCoords(NSolBackwards(j),k,2)
                zjb = ZSolCoords(NSolBackwards(j),k,3)
              else 
                xjf = xj
                yjf = yj
                zjf = zj
                xjb = xj
                yjb = yj
                zjb = zj
              endif
              if (.not. jperturbedatom) then
                dij2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c               Skip atom if no overlap
                if (dij2 .gt. (ri+rj)**2) goto 40
c               Eliminate buried elements
                do 50 l=1,nsurfpoints
c                 skip element if already buried (by another atom then)
                  if (iskip(l)) goto 50
                  dkj2 = ( ipoint(3*l-2)-xj )**2 + ( ipoint(3*l-1)-yj )**2 +
     &              ( ipoint(3*l) - zj)**2
                  if (dkj2 .lt. rj2) then 
                    iskip(l) = .true.
                    iskipf(l) = .true.
                    iskipb(l) = .true.
                  endif
 50             continue            
              else
c               J is a fep atom and must do fwd/bwd/ref independently
                dij2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
                if (jvanish .or. dij2 .gt. (ri+rj)**2) goto 510
                do 51 l=1,nsurfpoints
                  if (iskip(l)) goto 51
                  dkj2 = ( ipoint(3*l-2)-xj )**2 + ( ipoint(3*l-1)-yj )**2 +
     &              ( ipoint(3*l) - zj)**2
                  if (dkj2 .lt. rj2) iskip(l) = .true.                
 51             continue
 510            continue
c               Fwd
                dijf2 = (xif-xjf)**2 + (yif-yjf)**2 + (zif-zjf)**2
                if (jvanishF .or. dijf2 .gt. (rif+rjf)**2) goto 520
                do 52 l=1,nsurfpoints
                  if (iskipf(l)) goto 52
                  dkj2 = ( ipointf(3*l-2)-xjf )**2 + ( ipointf(3*l-1)-yjf )**2 +
     &              ( ipointf(3*l) - zjf)**2
                  if (dkj2 .lt. rjf2) iskipf(l) = .true.
 52             continue
 520            continue                
c               Bwd
                dijb2 = (xib-xjb)**2 + (yib-yjb)**2 + (zib-zjb)**2
                if (jvanishB .or. dijb2 .gt. (rib+rjb)**2) goto 530
                do 53 l=1,nsurfpoints
                  if (iskipb(l)) goto 53
                  dkj2 = ( ipointb(3*l-2)-xjb )**2 + ( ipointb(3*l-1)-yjb )**2 +
     &              ( ipointb(3*l) - zjb)**2
                  if (dkj2 .lt. rjb2) iskipb(l) = .true.
 53             continue
 530            continue
              endif
 40         continue              
 30       continue
        endif
c       Sum remaining elements
        isurf = 0.0d0
        isurff = 0.0d0
        isurfb = 0.0d0
        do 200 k=1,nsurfpoints
          if (.not. iskip(k) ) isurf = isurf + darea
          if (.not. iskipf(k) ) isurff = isurff + dareaf
          if (.not. iskipb(k) ) isurfb = isurfb + dareab
 200    continue
c       write(printstring,*) 'svn ',isvn,' at ',i,' nel ',isurf/darea
c       call printLine(DEBUG,printstring)
C       Increment the total surface area
        oldsurf = OldSvnAtomSurf(isvn,i)
        oldsurff = OldSvnAtomSurfF(isvn,i)
        oldsurfb = OldSvnAtomSurfB(isvn,i)

        diffsurf = isurf - oldsurf
        diffsurfF  = isurff - oldsurff
        diffsurfB = isurfb - oldsurfb        
C       Save isurf
        if (savenrg .eq. NEW) then
          NewSvnAtomSurf(isvn,i) = isurf
          NewSvnAtomSurfF(isvn,i) = isurff
          NewSvnAtomSurfB(isvn,i) = isurfb
        else if (savenrg .eq. OLD) then
          OldSvnAtomSurf(isvn,i) = isurf
          OldSvnAtomSurfF(isvn,i) = isurff
          OldSvnAtomSurfB(isvn,i) = isurfb
        endif

        SolventSurface = SolventSurface + diffsurf
        sanrg = sanrg + diffsurf * isurftens
        if (DifferentLambdas) then
          surfacef = surfacef + diffsurfF
          surfaceb = surfaceb + diffsurfB
          sanrgf = sanrgf + diffsurfF * isurftensf
          sanrgb = sanrgb + diffsurfB * isurftensb
        else
          surfacef = SolventSurface
          surfaceb = SolventSurface
          sanrgf = sanrg
          sanrgb = sanrg
        endif 
c        write (printstring,85) isvn,i,isurf,isurff,isurfb
c        call printLine(DEBUG,printstring) 
c 85     format('Solvent ',i3,' Atom ',i4,' SA ',f7.3,
c     &    ' SAF ',f7.3, ' SAB ',f7.3)
c        write (printstring,86) isvn,i,diffsurf*isurftens,diffsurfF * isurftensf,
c     .    diffsurfB * isurftensb
c 86     format('Solvent ',i3,' Atom ',i4,' SAnrg ',f7.4,
c     &    ' SAnrgF ',f7.4, ' SAnrgB ',f7.4)
c        call printLine(DEBUG,printstring) 
c       Add the energy 
c        if (SvnAtSAcount(isvn,i)) fastSolventSurfaceShrake = 
c     &    fastSolventSurfaceShrake + diffsurf * SURFTENS

 10   continue
c$$$      NewSolventSurf(isvn) = OldSolventSurf(isvn) + SolventSurface
c$$$      NewSolventSurfF(isvn) = OldSolventSurfF(isvn) + surfacef
c$$$      NewSolventSurfB(isvn) = OldSolventSurfB(isvn) + surfaceb
c      write(printstring,86) isvn,SolventSurface,surfacef,surfaceb
c      call printLine(DEBUG,printstring)
c 86   format('Solvent ',i3,' Surface Area ',f8.3, ' Fwd ',f8.3,' Bwd ',f8.3)
      return
      end
