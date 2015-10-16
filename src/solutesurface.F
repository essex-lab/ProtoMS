c
c     To avoid bugs in the calculation of the Surface Area, it is important
c     that the XX--Solute is done last to make sure that, for a non fep atom,
c     if iskip(m) is set to true, then iskipf(m) and iskipb(m) were previously 
c     set to true.
c     God I hope I will still understand this in 3 months

      double precision function SoluteSurface(savenrg,isol,surfacef,surfaceb,sanrg,sanrgf,sanrgb)
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
      integer isol
      integer i,j,k,l,m,nats,natsother,nres,resstart,resend
      double precision ipoint(nsurfpoints*3)
      double precision iarea(nsurfpoints),darea
      logical  iskip(nsurfpoints)
      double precision ri,ri2,rj,rj2,dij2,dkj2
      double precision xi,yi,zi,xj,yj,zj
      double precision isurf,isurftens
      integer ipar,jpar
      double precision rif,rif2,isurftensf,rib,rib2,isurftensb
      double precision sanrg,sanrgf,sanrgb
      double precision rjf,rjf2,rjb,rjb2,surfacef,surfaceb,dijf2,dijb2
      double precision xif,yif,zif,xib,yib,zib
      double precision xjf,yjf,zjf,xjb,yjb,zjb
      logical iperturbedatom,jperturbedatom
      double precision ipointf(nsurfpoints*3),iareaf(nsurfpoints),dareaf
      logical iskipf(nsurfpoints)
      double precision ipointb(nsurfpoints*3),iareab(nsurfpoints),dareab
      logical iskipb(nsurfpoints)
      double precision isurff,isurfb
      logical ivanish,ivanishF,ivanishB
      logical jvanish,jvanishF,jvanishB
      double precision oldsurf,diffsurf
      double precision oldsurff,oldsurfb,diffsurfF,diffsurfB
      double precision scl,sclf,sclb,sclj,sclbj,sclfj
      double precision getLowScale,getHighScale
      logical isolDual,jsoldual
      integer temp,tempf,tempb

      SoluteSurface = ZERO
      surfacef = ZERO
      surfaceb = ZERO
      sanrg = ZERO
      sanrgf = ZERO
      sanrgb = ZERO

      nats = NSolAtoms(isol)
c     see if the solute is part of a dual topology solute
      if (DualTopologySolute(isol).gt.0) then
c       this is the low dual topology solute
        isolDual = .true.
        scl = getLowScale(Lambda)
        sclb = getLowScale(LambdaB)
        sclf = getLowScale(LambdaF)
      else if (DualTopologySolute(isol).eq.0) then
c       this is not a dual topology solute
        isoldual = .false.
        scl = ONE
        sclf = ONE
        sclb = ONE
      else
c       this is the high dual topology solute
        isoldual = .true.
        scl = getHighScale(Lambda)
        sclb = getHighScale(LambdaB)
        sclf = getHighScale(LambdaF)
      endif

      do 10 i=1,nats
c       For each atom i
c       Check if it is necessary to compute the surface again or if it can be
c       regenerated
        if (.not. SolAtredoSA(isol,i)) then
c          isurf = OldSolAtomSurf(isol,i)
c          write(printstring,*) ' sol ',isol,' atom ',i,
c     &      'SA did not change, regenerated ',isurf
c          call printLine(DEBUG,printstring)
c          goto 201
          goto 10
        endif
c       Otherwise we have to do it again...
c       Skip atom if radii is 0
        iperturbedatom = .false.
        ivanish = .false.
        ivanishF = .false.
        ivanishB = .false.
        ipar = ZSolPar(isol,i)
        if (ipar .lt. 0) then
          iperturbedatom = .true.
c         Perturbed Atom
          ri = VariableSurfRadius(-ipar,2)
          isurftens = VariableSurfTension(-ipar,2)
          if (DifferentLambdas) then
c           Set up forward/backward ri and surface tensions
            rif = VariableSurfRadius(-ipar,5)
            isurftensf = VariableSurfTension(-ipar,5)
            rib = VariableSurfRadius(-ipar,4)
            isurftensb = VariableSurfTension(-ipar,4)
          else
            rif = ri
            isurftensf = isurftens
            rib = ri 
            isurftensb = isurftens
          endif
        else
c         Normal Atom
          ri = SurfRadDB(CLJType(ipar))
          isurftens = SurfTensionDB(CLJType(ipar))
          rif = ri
          isurftensf = isurftens
          rib = ri
          isurftensb = isurftens
        endif
c          write(printstring,*) 'ipar ',ipar,' ri ',ri,' rif ',rif,
c     &    ' rib ',rib,' isurftens ',isurftens, 'sftF ',isurftensf,
c     &    ' isurftensb ',isurftensb
c          call printLine(DEBUG,printstring)
c       Check if in one state the atom has vanished
c        if (iperturbedatom) then
        if (ri .lt. 0.01d0 .or. scl .lt. TINY) ivanish = .true.
        if (rif .lt. 0.01d0 .or. sclf .lt. TINY) ivanishF = .true.
        if (rib .lt. 0.01d0 .or. sclb .lt. TINY) ivanishB = .true.
        
c         
c If atom does not exist in any state then skip it completely
c
          if (ivanish .and. ivanishF .and. ivanishB) goto 10
c        else 
c         Normal atom with 0 radii (geometrical Dummy?)
c          if (ri .lt. 0.01d0) goto 10
c        endif
        ri = ri + probe
        ri2 = ri **2
        xi = ZSolCoords(isol,i,1)
        yi = ZSolCoords(isol,i,2)
        zi = ZSolCoords(isol,i,3)
        darea = ( FOURPI * ri2 ) / nsurfpoints
c        if (iperturbedatom .or. DualTopologySolute(isol) .ne. 0) then
        if (NSolVariable(isol) .gt. 0) then
c         FEP atom and coordinates lambda dependent
c         forward coordinates
          iperturbedatom = .true.
          xif = ZSolCoords(NSolForwards(isol),i,1)
          yif = ZSolCoords(NSolForwards(isol),i,2)
          zif = ZSolCoords(NSolForwards(isol),i,3)
c         backward coordinates
          xib = ZSolCoords(NSolBackwards(isol),i,1)
          yib = ZSolCoords(NSolBackwards(isol),i,2)
          zib = ZSolCoords(NSolBackwards(isol),i,3)
        else
c         FEP atom, coordinates not coupled to lambda
          xif = xi
          yif = yi
          zif = zi
          xib = xi
          yib = yi
          zib = zi
        endif
        rif = rif + probe
        rif2 = rif ** 2
        dareaf = ( FOURPI * rif2 ) / nsurfpoints
        rib = rib + probe
        rib2 = rib ** 2
        dareab = ( FOURPI * rib2 ) / nsurfpoints
c       else
c       Normal atom, forward==reference==backward
c       xif = xi
c       yif = yi
c       zif = zi
c       rif = ri
c       rif2 = ri2
c       dareaf = darea
c       xib = xi
c       yib = yi
c       zib = zi
c       rib = ri
c       rib2 = ri2
c       dareab = darea
c       endif
c       write(printstring,*) 'xi ',xi,' yi ',yi,' zi ',zi,' ri ',ri
c       call printLine(DEBUG,printstring)
        
        if (isoldual) iperturbedatom = .true.

c       Set up the surface elements of atom i
        if (SAquality .eq. 1) then
          do 20 j=1,nsurfpoints
            ipoint(3*j-2) = points12(3*j-2) * ri + xi
            ipoint(3*j-1) = points12(3*j-1) * ri + yi
            ipoint(3*j) = points12(3*j) * ri + zi
            iarea(j) = darea
            iskip(j) = .false.
 20       continue
c         Set up forward/backward surface elements
          if (iperturbedatom) then
            do 201 j=1,nsurfpoints
c             forward
              ipointf(3*j-2) = points12(3*j-2) * rif + xif
              ipointf(3*j-1) = points12(3*j-1) * rif + yif
              ipointf(3*j) = points12(3*j) * rif + zif
              iareaf(j) = dareaf
              iskipf(j) = .false.
c             backward
              ipointb(3*j-2) = points12(3*j-2) * rib + xib
              ipointb(3*j-1) = points12(3*j-1) * rib + yib
              ipointb(3*j) = points12(3*j) * rib + zib
              iareab(j) = dareab
              iskipb(j) = .false.
 201        continue
          else
            do 202 j=1,nsurfpoints
c             forward
              ipointf(3*j-2) = ipoint(3*j-2) 
              ipointf(3*j-1) = ipoint(3*j-1) 
              ipointf(3*j) = ipoint(3*j) 
              iareaf(j) = dareaf
              iskipf(j) = .false.
c             backward
              ipointb(3*j-2) = ipoint(3*j-2) 
              ipointb(3*j-1) = ipoint(3*j-1) 
              ipointb(3*j) = ipoint(3*j) 
              iareab(j) = dareab
              iskipb(j) = .false.
 202        continue
          endif              
        else if (SAquality .eq. 2) then
          do 21 j=1,nsurfpoints
            ipoint(3*j-2) = points42(3*j-2) * ri + xi
            ipoint(3*j-1) = points42(3*j-1) * ri + yi
            ipoint(3*j) = points42(3*j) * ri + zi
            iarea(j) = darea
            iskip(j) = .false.
 21       continue    
          if (iperturbedatom) then
            do 211 j=1,nsurfpoints
c             forward
              ipointf(3*j-2) = points42(3*j-2) * rif + xif
              ipointf(3*j-1) = points42(3*j-1) * rif + yif
              ipointf(3*j) = points42(3*j) * rif + zif
              iareaf(j) = dareaf
              iskipf(j) = .false.
c             backward
              ipointb(3*j-2) = points42(3*j-2) * rib + xib
              ipointb(3*j-1) = points42(3*j-1) * rib + yib
              ipointb(3*j) = points42(3*j) * rib + zib
              iareab(j) = dareab
              iskipb(j) = .false.
 211        continue
          else
            do 212 j=1,nsurfpoints
c             forward
              ipointf(3*j-2) = ipoint(3*j-2) 
              ipointf(3*j-1) = ipoint(3*j-1) 
              ipointf(3*j) = ipoint(3*j) 
              iareaf(j) = dareaf
              iskipf(j) = .false.
c             backward
              ipointb(3*j-2) = ipoint(3*j-2) 
              ipointb(3*j-1) = ipoint(3*j-1) 
              ipointb(3*j) = ipoint(3*j) 
              iareab(j) = dareab
              iskipb(j) = .false.
 212        continue
          endif      
        else if (SAquality .eq. 3) then
          do 22 j=1,nsurfpoints
            ipoint(3*j-2) = points162(3*j-2) * ri + xi
            ipoint(3*j-1) = points162(3*j-1) * ri + yi
            ipoint(3*j) = points162(3*j) * ri + zi
            iarea(j) = darea
            iskip(j) = .false.
 22       continue 
          if (iperturbedatom) then
            do 222 j=1,nsurfpoints
c             forward
              ipointf(3*j-2) = points162(3*j-2) * rif + xif
              ipointf(3*j-1) = points162(3*j-1) * rif + yif
              ipointf(3*j) = points162(3*j) * rif + zif
              iareaf(j) = dareaf
              iskipf(j) = .false.
c             backward
              ipointb(3*j-2) = points162(3*j-2) * rib + xib
              ipointb(3*j-1) = points162(3*j-1) * rib + yib
              ipointb(3*j) = points162(3*j) * rib + zib
              iareab(j) = dareab
              iskipb(j) = .false.
 222        continue
          else
            do 223 j=1,nsurfpoints
c             forward
              ipointf(3*j-2) = ipoint(3*j-2) 
              ipointf(3*j-1) = ipoint(3*j-1) 
              ipointf(3*j) = ipoint(3*j) 
              iareaf(j) = dareaf
              iskipf(j) = .false.
c             backward
              ipointb(3*j-2) = ipoint(3*j-2) 
              ipointb(3*j-1) = ipoint(3*j-1) 
              ipointb(3*j) = ipoint(3*j) 
              iareab(j) = dareab
              iskipb(j) = .false.
 223        continue
          endif         
        else if (SAquality .eq. 4) then
          do 23 j=1,nsurfpoints
            ipoint(3*j-2) = points642(3*j-2) * ri + xi
            ipoint(3*j-1) = points642(3*j-1) * ri + yi
            ipoint(3*j) = points642(3*j) * ri + zi
            iarea(j) = darea
            iskip(j) = .false.
 23       continue     
          if (iperturbedatom) then
            do 233 j=1,nsurfpoints
c             forward
              ipointf(3*j-2) = points642(3*j-2) * rif + xif
              ipointf(3*j-1) = points642(3*j-1) * rif + yif
              ipointf(3*j) = points642(3*j) * rif + zif
              iareaf(j) = dareaf
              iskipf(j) = .false.
c             backward
              ipointb(3*j-2) = points642(3*j-2) * rib + xib
              ipointb(3*j-1) = points642(3*j-1) * rib + yib
              ipointb(3*j) = points642(3*j) * rib + zib
              iareab(j) = dareab
              iskipb(j) = .false.
 233        continue
          else
            do 234 j=1,nsurfpoints
c             forward
              ipointf(3*j-2) = ipoint(3*j-2) 
              ipointf(3*j-1) = ipoint(3*j-1) 
              ipointf(3*j) = ipoint(3*j) 
              iareaf(j) = dareaf
              iskipf(j) = .false.
c             backward
              ipointb(3*j-2) = ipoint(3*j-2) 
              ipointb(3*j-1) = ipoint(3*j-1) 
              ipointb(3*j) = ipoint(3*j) 
              iareab(j) = dareab
              iskipb(j) = .false.
 234        continue
          endif     
        endif
c       Now do Solute with other Proteins
        if (NProteins .gt. ZERO) then
          do 60 j=1,NProteins
            nres = NResidues(j)
            do 70 k=1,nres
              resstart = BBZid(j,k,1)
              resend = BBZid(j,k,4)
              do 80 l=resstart,resend
                jpar = ZPar(j,l)
                if (jpar .lt. 0) then
                  rj = 0.0d0
                  rj = rj + probe
                  rj2 = rj**2
                  rjf = rj
                  rjb = rj
                  rjf2 = rj2
                  rjb2 = rj2
                  write(printstring,*) 'Cannot do Shrake of perturbed ',
     &              ' protein atom  ' 
                  call printLine(FATAL,printstring)
                  call closeProgram(-1)
                else
                  rj = SurfRadDB(CLJType(jpar))
c                 J non FEP, no SA (geometrical dummy?)
                  if (rj .lt. 0.01d0) goto 80
                  rj = rj + probe
                  rj2 = rj**2
                  rjf = rj
                  rjb = rj
                  rjf2 = rj2
                  rjb2 = rj2
                endif
                xj = ZCoords(j,l,1)
                yj = ZCoords(j,l,2)
                zj = ZCoords(j,l,3)
c               I don't think you can do FEP on Protein at the moment !
                xjf = xj
                yjf = yj
                zjf = zj
                xjb = xj
                yjb = yj
                zjb = zj
                if (.not. iperturbedatom) then
                  dij2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c                 Skip atom if no overlap
                  if (dij2 .gt. (ri+rj)**2) goto 80
c                 Eliminate buried elements
                  do 90 m=1,nsurfpoints
c                   skip element if already buried (by another atom then)
                    if (iskip(m)) goto 90
                    dkj2 = ( ipoint(3*m-2)-xj )**2 + ( ipoint(3*m-1)-yj )**2 +
     &                ( ipoint(3*m) - zj)**2
                    if (dkj2 .lt. rj2) then 
                      iskip(m) = .true.
                      iskipf(m) = .true.
                      iskipb(m) = .true.
                    endif
 90               continue
                else
c                 Doing fwd/bwd because the solute atom is FEP
c                 Do the reference calculation if there is some overlap
                  dij2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c                 Skip if one of the atoms in this state has 'vanished' or if there 
c                 is no overlap
                  if (ivanish .or. dij2 .gt. (ri+rj)**2) goto 610
                  do 91 m=1,nsurfpoints
                    if (iskip(m)) goto 91
                    dkj2 = ( ipoint(3*m-2)-xj )**2 + ( ipoint(3*m-1)-yj )**2 +
     &                ( ipoint(3*m) - zj)**2
                    if (dkj2 .lt. rj2) then 
                      iskip(m) = .true.
                    endif
 91               continue
 610              continue
c                 Do the same job with forward
                  dijf2 = (xif-xjf)**2 + (yif-yjf)**2 + (zif-zjf)**2
                  if (ivanishF .or. dijf2 .gt. (rif+rjf)**2) goto 620
                  do 92 m=1,nsurfpoints
                    if (iskipf(m)) goto 92
                    dkj2 = ( ipointf(3*m-2)-xjf )**2 + ( ipointf(3*m-1)-yjf )**2 +
     &                ( ipointf(3*m) - zjf)**2
                    if (dkj2 .lt. rjf2) then 
                      iskipf(m) = .true.
                    endif
 92               continue
 620              continue
c                 Now backward
                  dijb2 = (xib-xjb)**2 + (yib-yjb)**2 + (zib-zjb)**2
                  if (ivanishB .or. dijb2 .gt. (rib+rjb)**2) goto 630
                  do 93 m=1,nsurfpoints
                    if (iskipb(m)) goto 93
                    dkj2 = ( ipointb(3*m-2)-xjb )**2 + ( ipointb(3*m-1)-yjb )**2 +
     &                ( ipointb(3*m) - zjb)**2
                    if (dkj2 .lt. rjb2) then 
                      iskipb(m) = .true.
                    endif
 93               continue
 630              continue                  
                endif
 80           continue
 70         continue
 60       continue
        endif
c       Now do Solute with other solvents
        if (NSolvents .gt. ZERO) then
          do 100 j=1,NSolvents
            natsother = NSvnAtoms(j)
            do 110 k=1,natsother
              jpar = ZSvnPar(j,k)
              if (jpar .lt. 0) then
                rj = 0.0d0
                rj = rj + probe
                rj2 = rj**2
                rjf = rj
                rjb = rj
                rjf2 = rj2
                rjb2 = rj2
                write(printstring,*) 'Cannot do Shrake of perturbed ',
     &            ' solvent atom  ' 
                call printLine(FATAL,printstring)
                call closeProgram(-1)
              else
                rj = SurfRadDB(CLJType(jpar))
c               J non FEP, no SA (geometrical dummy?)
                if (rj .lt. 0.01d0) goto 110
                rj = rj + probe
                rj2 = rj**2
                rjf = rj
                rjb = rj
                rjf2 = rj2
                rjb2 = rj2
              endif              
              xj = ZSvnCoords(j,k,1)
              yj = ZSvnCoords(j,k,2)
              zj = ZSvnCoords(j,k,3)
              xjf = xj
              yjf = yj
              zjf = zj
              xjb = xj
              yjb = yj
              zjb = zj
              if (.not. iperturbedatom) then
                dij2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c               Skip atom if no overlap
                if (dij2 .gt. (ri+rj)**2) goto 110
c               Eliminate buried elements
                do 120 l=1,nsurfpoints
c                 skip element if already buried (by another atom then)
                  if (iskip(l)) goto 120
                  dkj2 = ( ipoint(3*l-2)-xj )**2 + ( ipoint(3*l-1)-yj )**2 +
     &              ( ipoint(3*l) - zj)**2
                  if (dkj2 .lt. rj2) then 
                    iskip(l) = .true.
                    iskipf(l) = .true.
                    iskipb(l) = .true.
                  endif
 120            continue
              else
c               Doing fwd/bwd because the solute atom is FEP
c               Do the reference calculation if there is some overlap
                dij2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c               Skip if one of the atoms in this state has 'vanished' or if there 
c               is no overlap
                if (ivanish .or. dij2 .gt. (ri+rj)**2) goto 710
                do 121 l=1,nsurfpoints
                  if (iskip(l)) goto 121
                  dkj2 = ( ipoint(3*l-2)-xj )**2 + ( ipoint(3*l-1)-yj )**2 +
     &              ( ipoint(3*l) - zj)**2
                  if (dkj2 .lt. rj2) then 
                    iskip(l) = .true.
                  endif
 121            continue
 710            continue
c               Do the same job with forward
                dijf2 = (xif-xjf)**2 + (yif-yjf)**2 + (zif-zjf)**2
                if (ivanishF .or. dijf2 .gt. (rif+rjf)**2) goto 720
                do 122 l=1,nsurfpoints
                  if (iskipf(l)) goto 122
                  dkj2 = ( ipointf(3*l-2)-xjf )**2 + ( ipointf(3*l-1)-yjf )**2 +
     &              ( ipointf(3*l) - zjf)**2
                  if (dkj2 .lt. rjf2) then 
                    iskipf(l) = .true.
                  endif
 122            continue
 720            continue
c               Now backward
                dijb2 = (xib-xjb)**2 + (yib-yjb)**2 + (zib-zjb)**2
                if (ivanishB .or. dijb2 .gt. (rib+rjb)**2) goto 730
                do 123 l=1,nsurfpoints
                  if (iskipb(l)) goto 123
                  dkj2 = ( ipointb(3*l-2)-xjb )**2 + ( ipointb(3*l-1)-yjb )**2 +
     &              ( ipointb(3*l) - zjb)**2
                  if (dkj2 .lt. rjb2) then 
                    iskipb(l) = .true.
                  endif
 123            continue
 730            continue                  
              endif
 110        continue
 100      continue
        endif
c       Now do Solute with other Solutes
        do 30 j=1,NSolutes
          if (IgnoreSolute(j)) goto 30
c         If it is a pair of dual topology solute then they should ignore each other ! 
          if (DualTopologySolute(isol).eq. j .or. DualTopologySolute(j).eq. isol) goto 30
c         Get the scaling factors of that solute, to decide if some atoms must be skipped
          if (DualTopologySolute(j).gt.0) then
c           this is the low dual topology solute
            jsolDual = .true.
            sclj = getLowScale(Lambda)
            sclbj = getLowScale(LambdaB)
            sclfj = getLowScale(LambdaF)
          else if (DualTopologySolute(j).eq.0) then
c           this is not a dual topology solute
            jsoldual = .false.
            sclj = ONE
            sclfj = ONE
            sclbj = ONE
          else
c           this is the high dual topology solute
            jsoldual = .true.
            sclj = getHighScale(Lambda)
            sclbj = getHighScale(LambdaB)
            sclfj = getHighScale(LambdaF)
          endif          
          natsother = NSolAtoms(j)
          do 40 k=1,natsother
C           Skip itself !
            if (i .eq. k .and. isol. eq. j) goto 40
            jperturbedatom = .false.
            jvanish = .false.
            jvanishF = .false.
            jvanishB = .false.
            jpar = ZSolPar(j,k)
            if (jpar .lt. 0) then
              jperturbedatom = .true.
c             Perturbed Atom
              rj = VariableSurfRadius(-jpar,2)
              if (DifferentLambdas) then
c               Set up forward/backward ri and surface tensions
                rjf = VariableSurfRadius(-jpar,5)
                rjb = VariableSurfRadius(-jpar,4)
              else
                rjf = rj
                rjb = rj 
              endif
            else
c             Normal Atom
              rj = SurfRadDB(CLJType(jpar))
              rjf = rj
              rjb = rj
            endif

            if (rj  .lt. 0.01d0 .or. sclj  .lt. TINY ) jvanish = .true.
            if (rjf .lt. 0.01d0 .or. sclfj .lt. TINY ) jvanishF = .true.
            if (rjb .lt. 0.01d0 .or. sclbj .lt. TINY ) jvanishB = .true.
            rj = rj + probe
            rj2 = rj **2
            rjf = rjf + probe
            rjf2 = rjf ** 2
            rjb = rjb + probe
            rjb2 = rjb ** 2
            if (jsoldual) jperturbedatom = .true.
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
c            write(printstring,*) 'isol ',isol,' i ',i,' with jsol ',j,
c     &        ' j ',k,' ivan ',ivanish,ivanishF,ivanishB,
c     &        ' jvan ',jvanish,jvanishF,jvanishB
c            call printLine(DEBUG,printstring)
c            write(printstring,*) 'ipert ',iperturbedatom,' jpert ',jperturbedatom
c            call printLine(DEBUG,printstring)

            if (.not. iperturbedatom .and. .not. jperturbedatom) then
              dij2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c             Skip atom if no overlap
              if (ivanish .or. jvanish .or. dij2 .gt. (ri+rj)**2 ) goto 40
c             Eliminate buried elements
              do 50 l=1,nsurfpoints
c               skip element if already buried (by another atom then)
                if (iskip(l)) then
                  goto 50
                endif
                dkj2 = ( ipoint(3*l-2)-xj )**2 + ( ipoint(3*l-1)-yj )**2 +
     &            ( ipoint(3*l) - zj)**2
                if (dkj2 .lt. rj2) then 
                  iskip(l) = .true.
c                 Since both atoms are non FEP, we don't have to do the 
c                 forward/backward calculation to find this
                  iskipf(l) = .true.
                  iskipb(l) = .true.
                endif
 50           continue            
            else 
c             Otherwise..one or two atoms are FEP and we must do forward/backward
c             independently
c             Do the reference calculation if there is some overlap
              dij2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c             Skip if one of the atoms in this state has 'vanished' or if there 
c             is no overlap
              if (ivanish .or. jvanish .or. dij2 .gt. (ri+rj)**2) goto 510
              do 51 l=1,nsurfpoints
                if (iskip(l)) goto 51
                dkj2 = ( ipoint(3*l-2)-xj )**2 + ( ipoint(3*l-1)-yj )**2 +
     &            ( ipoint(3*l) - zj)**2
                if (dkj2 .lt. rj2) iskip(l) = .true.                
 51           continue
 510          continue
c             Do the same job with forward
              dijf2 = (xif-xjf)**2 + (yif-yjf)**2 + (zif-zjf)**2
              if (ivanishF .or. jvanishF .or. dijf2 .gt. (rif+rjf)**2) goto 520
              do 52 l=1,nsurfpoints
                if (iskipf(l)) goto 52
                dkj2 = ( ipointf(3*l-2)-xjf )**2 + ( ipointf(3*l-1)-yjf )**2 +
     &            ( ipointf(3*l) - zjf)**2
                if (dkj2 .lt. rjf2) iskipf(l) = .true.
 52           continue
 520          continue
c             Now backward
              dijb2 = (xib-xjb)**2 + (yib-yjb)**2 + (zib-zjb)**2
              if (ivanishB .or. jvanishB .or. dijb2 .gt. (rib+rjb)**2) goto 530
              do 53 l=1,nsurfpoints
                if (iskipb(l)) goto 53
                dkj2 = ( ipointb(3*l-2)-xjb )**2 + ( ipointb(3*l-1)-yjb )**2 +
     &            ( ipointb(3*l) - zjb)**2
                if (dkj2 .lt. rjb2) iskipb(l) = .true.
 53           continue
 530          continue
            endif
            temp = 0
            tempf = 0
            tempb = 0
            do 888 l=1,nsurfpoints
              if (.not. iskip(l)) temp = temp + 1
              if (.not. iskipf(l)) tempf = tempf + 1
              if (.not. iskipb(l)) tempb = tempb + 1
 888        continue
c            write(printstring,*) 'Non skipped ',temp,' F ',tempf,' B ',tempB
c            call printLine(DEBUG,printstring)
 40       continue
 30     continue
c       Sum remaining elements
        isurf = 0.0d0
        isurff = 0.0d0
        isurfb = 0.0d0
        do 200 k=1,nsurfpoints
          if (.not. iskip(k) )  isurf = isurf + darea
          if (.not. iskipf(k) ) isurff = isurff + dareaf
          if (.not. iskipb(k) ) isurfb = isurfb + dareab
 200    continue
c        Is zero if one atom has vanished
         if (ivanish) isurf = 0.0d0
         if (ivanishF) isurff = 0.0d0
         if (ivanishB) isurfb = 0.0d0        
C       Increment the total surface area
         oldsurf = OldSolAtomSurf(isol,i)
         oldsurff = OldSolAtomSurfF(isol,i)
         oldsurfb = OldSolAtomSurfB(isol,i)
         diffsurf = isurf - oldsurf
         diffsurfF  = isurff - oldsurff
         diffsurfB = isurfb - oldsurfb
C        Save isurf
         if (savenrg .eq. NEW) then
           NewSolAtomSurf(isol,i) = isurf
           NewSolAtomSurfF(isol,i) = isurff
           NewSolAtomSurfB(isol,i) = isurfb
         else if (savenrg .eq. OLD) then
c          This should only happen during initialisation...
           OldSolAtomSurf(isol,i) = isurf
           OldSolAtomSurfF(isol,i) = isurff
           OldSolAtomSurfB(isol,i) = isurfb
         endif

c       Add the new - old energy 
c        if (SolAtSAcount(isol,i)) fastSoluteSurfaceShrake = 
c     &    fastSoluteSurfaceShrake + diffsurf * SURFTENS
c        write(printstring,*) 'OldSurf ',oldsurf,' NewSurf ',isurf,
c     &    ' diffsurf ',diffsurf,' isurftens ',isurftens
c        call printLine(DEBUG,printstring)
c       write(printstring,*) ' tens ',isurftens,' tensF ',isurftensF,' tensB ',isurftensB
c        call printLine(DEBUG,printstring)
        SoluteSurface = SoluteSurface + diffsurf
        sanrg = sanrg + diffsurf * isurftens
        if (DifferentLambdas) then
          surfacef = surfacef + diffsurfF
          surfaceb = surfaceb + diffsurfB
          sanrgf = sanrgf + diffsurfF * isurftensf
          sanrgb = sanrgb + diffsurfB * isurftensb
        else
          surfacef = SoluteSurface
          surfaceb = SoluteSurface
          sanrgf = sanrg
          sanrgb = sanrg
        endif
c        write (printstring,85) isol,i,isurf,isurff,isurfb
c        call printLine(DEBUG,printstring) 
c 85     format('Solute ',i3,' Atom ',i4,' SA ',f8.3,' SAF ',f8.3,' SAB ',f8.3)
 10   continue

c      NewSoluteSurf(isol) = OldSoluteSurf(isol) + SoluteSurface
c      NewSoluteSurfF(isol) = OldSoluteSurfF(isol) + surfacef
c      NewSoluteSurfB(isol) = OldSoluteSurfB(isol) + surfaceb
c      write(printstring,86) isol,SoluteSurface,surfacef,surfaceb
c      call printLine(DEBUG,printstring)
c 86   format('Solute ',i3,' Surface Area ',f8.3,' Fwd ',f8.3,' Bwd ',f8.3)

      return
      end
