
      subroutine markSAupdate
      include 'dimensions.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'solvents.inc'
      include 'energies.inc'
      include 'flags.inc'
      include 'constants.inc'
      include 'printer.inc'
      include 'surface.inc'      
c###############################################
c     Determines which part of the system needs
c     to have its Surface Area updated
c
c     Julien Michel November 2004
c###############################################
      integer i,j,k,l,m,m2
      logical hasSoluteMoved
      logical hasResidueMoved
      logical hasSolventMoved
      integer nats,nats2,nres,resstart,resend,nres2,resstart2,resend2
      integer movedsol,movedsvn,movedres
      double precision xi,yi,zi,dij2
      double precision xif,yif,zif,dijf2
      double precision xib,yib,zib,dijb2
      logical dogeom
      double precision xj,yj,zj
      double precision xjf,yjf,zjf
      double precision xjb,yjb,zjb
c
c     Note : At the moment if any atom may need to be updated in a residue,
c     then every atom of the residue is done again. That can be improved
c     by doing only the atoms in the residue that needs updating...
c

c
c     First clear redoSA flags
c
      if (NSolutes .gt. 0) then
        do 1 i=1,NSolutes
          if (IgnoreSolute(i)) goto 1
          SoluteredoSA(i) = .false.
          nats = NSoLAtoms(i)
          do 2 j=i,nats
            SolAtredoSA(i,j) = .false.
 2        continue
 1     continue
      endif

      if (NSolvents .gt. 0) then
        do 3 i=1,NSolvents
          SolventredoSA(i) = .false.
          nats = NSvnAtoms(i)
          do 4 j=1,nats
            SvnAtredoSA(i,j) = .false.            
 4        continue
 3      continue
      endif

      if (NProteins .gt. 0) then
        do 5 i=1,NProteins
          nres = NResidues(i)
          do 6 j=1,nres
            ResidueredoSA(i,j) = .false.
            resstart = BBZid(i,j,1)
            resend = BBZid(i,j,4)            
            do 7 k=resstart,resend
              ProtAtredoSA(i,k) = .false.
 7          continue
 6        continue
 5      continue
      endif

c
c     Now set the flags
c

      if (NSolutes .gt. 0) then
        do 10 i=1,NSolutes
          if (IgnoreSolute(i)) goto 10
          if (hasSoluteMoved(i)) then
            movedsol = i
            SoluteredoSA(i) = .true.
            dogeom = .false.
            if (DifferentLambdas .and. NSolVariable(movedsol) .gt. 0) dogeom = .true.
c           Do pair pair of every atom in Solute with the rest of the system, mark
c           atoms which needs update
            nats = NSolAtoms(movedsol)
            do 20 j=1,nats
              SolAtredoSA(i,j) = .true.
              xi = ZSolCoords(movedsol,j,1)
              yi = ZSolCoords(movedsol,j,2)
              zi = ZSolCoords(movedsol,j,3)
              if (dogeom) then
                xif = ZSolCoords(NSolForwards(movedsol),j,1)
                yif = ZSolCoords(NSolForwards(movedsol),j,2)
                zif = ZSolCoords(NSolForwards(movedsol),j,3)
                xib = ZSolCoords(NSolBackwards(movedsol),j,1)
                yib = ZSolCoords(NSolBackwards(movedsol),j,2)
                zib = ZSolCoords(NSolBackwards(movedsol),j,3)
              else
                xif = xi
                yif = yi
                zif = zi
                xib = xi
                yib = yi
                zib = zi
              endif
c             Do Solute--Solute 
              do 30 k=1,NSolutes
                if (IgnoreSolute(k)) goto 30
                if (k .eq. movedsol) goto 30
c               If already marked for update just skip it
c               if (SoluteredoSA(k) .or. hasSoluteMoved(k)) goto 30
                nats2 = NSolAtoms(k)
                do 40 l=1,nats2
                  xj = ZSolCoords(k,l,1)
                  yj = ZSolCoords(k,l,2)
                  zj = ZSolCoords(k,l,3)
                  dij2 = (xi-xj)**2 + (yi-yj)**2 +(zi-zj)**2 
                  if (dogeom) then
c                   Check if solute k is also dogeom
                    if (NSolVariable(k) .gt. 0) then
c                     Solute k coords change with lambda
                      xjf = ZSolCoords(NSolForwards(k),l,1)
                      yjf = ZSolCoords(NSolForwards(k),l,2)
                      zjf = ZSolCoords(NSolForwards(k),l,3)
                      xjb = ZSolCoords(NSolBackwards(k),l,1)
                      yjb = ZSolCoords(NSolBackwards(k),l,2)
                      zjb = ZSolCoords(NSolBackwards(k),l,3)
                    else
c                     Solute k coords do not change with Lamdba
                      xjf = xj
                      yjf = yj
                      zjf = zj
                      xjb = xj
                      yjb = yj
                      zjb = zj
                    endif
c                   Do forward
                    dijf2 = (xif-xjf)**2 + (yif-yjf)**2 +(zif-zjf)**2 
c                   Do backward
                    dijb2 = (xib-xjb)**2 + (yib-yjb)**2 +(zib-zjb)**2 
                  else
                    dijf2 = dij2
                    dijb2 = dij2
                  endif
                  if (dij2 .lt. SurfCut2 .or. dijf2 .lt. SurfCut2 
     &              .or. dijb2 .lt. SurfCut2) then 
c                   Mark that this solute needs to be redone and go to next step
                    SoluteredoSA(k) = .true.
                    SolAtredoSA(k,l) = .true.
                  endif                  
 40             continue
 30           continue
c             Do Solute--Proteins
              if (NProteins .gt. 0) then
                do 50 k=1,NProteins
                  nres = NResidues(k)
                  do 60 l=1,nres
c                   if (ResidueredoSA(k,l) .or. hasResidueMoved(k,l)) goto 60
                    resstart = BBZid(k,l,1)
                    resend = BBZid(k,l,4)
                    do 70 m=resstart,resend
                      dij2 = ( (xi-ZCoords(k,m,1))**2 + 
     &                  (yi-ZCoords(k,m,2))**2 +
     &                  (zi-ZCoords(k,m,3))**2 )
                      if (dij2 .lt. SurfCut2) then 
                        ResidueredoSA(k,l) = .true.
                        ProtAtredoSA(k,m) = .true.
c                       goto 60
                      endif
 70                 continue
 60               continue
 50             continue
              endif
c             Do Solute--Solvents
              if (NSolvents .gt. 0) then
                do 80 k=1,NSolvents
c                  if (SolventredoSA(k) .or. hasSolventMoved(k)) goto 80
c                 Assume no need to redo
                  nats2 = NSvnAtoms(k)
                  do 90 l=1,nats2
                    dij2 = ( (xi-ZSvnCoords(k,l,1))**2 + 
     &                (yi-ZSvnCoords(k,l,2))**2 +
     &                (zi-ZSvnCoords(k,l,3))**2 )
                    if (dij2 .lt. SurfCut2) then 
                      SolventredoSA(k) = .true.
                      SvnAtredoSA(k,l) = .true.
c                      goto 80
                    endif
 90               continue
 80             continue
              endif
 20         continue
          endif
 10     continue
      endif

      if (NProteins .gt. 0) then
        do 100 i=1,NProteins
          nres = NResidues(i)
          do 110 j=1,nres
            if (hasResidueMoved(i,j)) then
              movedres = j
              ResidueredoSA(i,movedres) = .true.
              resstart = BBZid(i,movedres,1)
              resend = BBZid(i,movedres,4)
              do 120 k=resstart,resend
                ProtAtredoSA(i,k) = .true.
                xi = ZCoords(i,k,1)
                yi = ZCoords(i,k,2)
                zi = ZCoords(i,k,3)  
                xif = xi
                yif = yi
                zif = zi
                xib = xi
                yib = yi
                zib = zi
c               Do Residue--Residue
                do 130 l=1,NProteins
                  nres2 = NResidues(l)
                  do 140 m=1,nres2
                    if (m .eq. movedres .and. i .eq. l) goto 140
c                    if (ResidueredoSA(l,m) .or. hasResidueMoved(l,m)) goto 140
                    resstart2 = BBZid(l,m,1)
                    resend2 = BBZid(l,m,4)
                    do 150 m2=resstart2,resend2
                      dij2 = ( (xi-ZCoords(l,m2,1))**2 + 
     &                  (yi-ZCoords(l,m2,2))**2 +
     &                  (zi-ZCoords(l,m2,3))**2 )
                      if (dij2 .lt. SurfCut2) then 
                        ResidueredoSA(l,m) = .true.
                        ProtAtredoSA(l,m2) = .true.
c                        goto 140
                      endif
 150                continue
 140              continue
 130            continue
c               Do Residue--Solute
                if (NSolutes .gt. 0) then
                  do 160 l=1,NSolutes
                    if (IgnoreSolute(l)) goto 160
c                    if (SoluteredoSA(l) .or. hasSoluteMoved(l)) goto 160
                    nats = NSolAtoms(l)
                    do 170 m=1,nats
                      xj = ZSolCoordS(l,m,1)
                      yj = ZSolCoordS(l,m,2)
                      zj = ZSolCoordS(l,m,3)
                      dij2 = ( (xi-xj)**2 + (yi-yj)**2 +(zi-zj)**2 )
                      if (NSolVariable(l) .gt. 0) then
                        xjf = ZSolCoords(NSolForwards(l),m,1)
                        yjf = ZSolCoords(NSolForwards(l),m,2)
                        zjf = ZSolCoords(NSolForwards(l),m,3)
                        xjb = ZSolCoords(NSolBackwards(l),m,1)
                        yjb = ZSolCoords(NSolBackwards(l),m,2)
                        zjb = ZSolCoords(NSolBackwards(l),m,3)
                        dijf2 = (xif-xjf)**2 + (yif-yjf)**2 +(zif-zjf)**2 
                        dijb2 = (xib-xjb)**2 + (yib-yjb)**2 +(zib-zjb)**2 
                      else
                        xjf = xj
                        yjf = yj
                        zjf = zj
                        xjb = xj
                        yjb = yj
                        zjb = zj
                        dijf2 = dij2
                        dijb2 = dij2
                      endif

                      if (dij2 .lt. SurfCut2 .or. dijf2 .lt. SurfCut2 
     &                  .or. dijb2 .lt. SurfCut2) then 
                        SoluteredoSA(l) = .true.
                        SolAtredoSA(l,m) = .true.
c                        goto 160
                      endif
 170                continue
 160              continue
                endif
c               Do Residue-Solvents
                if (NSolvents .gt. 0) then
                  do 180 l=1,NSolvents
c                    if (SolventredoSA(l) .or. hasSolventmoved(l)) goto 180
                    nats = NSvnAtoms(l)
                    do 190 m=1,nats
                      dij2 = ( (xi-ZSvnCoords(l,m,1))**2 + 
     &                  (yi-ZSvnCoords(l,m,2))**2 +
     &                  (zi-ZSvnCoords(l,m,3))**2 )
                      if (dij2 .lt. SurfCut2) then 
                        SolventredoSA(l) = .true.
                        SvnAtredoSA(l,m) = .true.
c                        goto 180
                      endif
 190                continue
 180              continue
                endif
 120          continue
            endif
 110      continue
 100    continue
      endif
   
      if (NSolvents .gt. 0) then
        do 200 i=1,NSolvents
          if (hasSolventMoved(i) ) then
            movedsvn = i
            SolventredoSA(i) = .true.
c            write(printstring,*) 'This svn moved ',movedsvn
c            call printLine(DEBUG,printstring)
c           Do pair pair of every atom in Solute with the rest of the system, mark
c           atoms which needs update
            nats = NSvnAtoms(movedsvn)
            do 210 j=1,nats
              SvnAtredoSA(i,j) = .true.
              xi = ZSvnCoords(movedsvn,j,1)
              yi = ZSvnCoords(movedsvn,j,2)
              zi = ZSvnCoords(movedsvn,j,3)
c             Do Solvent--Solvents
              do 220 k=1,NSolvents
                if (k .eq. movedsvn) goto 220
c                if (SolventredoSA(k) .or. hasSolventmoved(k)) goto 220
                nats2 = NSvnAtoms(k)
                do 230 l=1,nats2
                  dij2 = ( (xi-ZSvnCoords(k,l,1))**2 + 
     &              (yi-ZSvnCoords(k,l,2))**2 +
     &              (zi-ZSvnCoords(k,l,3))**2 )
                  if (dij2 .lt. SurfCut2) then 
                    SolventredoSA(k) = .true.
                    SvnAtredoSA(k,l) = .true.
c                    goto 220
                  endif
 230            continue
 220          continue

c             Do Solvent--Solutes
              if (NSolutes .gt. 0) then
                do 240 k=1,NSolutes
                  if (IgnoreSolute(k)) goto 240
c                  if (SoluteredoSA(k) .or. hasSolutemoved(k)) goto 240
                  nats2 = NSolAtoms(k)
                  do 250 l=1,nats2
                    dij2 = ( (xi-ZSolCoords(k,l,1))**2 + 
     &                (yi-ZSolCoords(k,l,2))**2 +
     &                (zi-ZSolCoords(k,l,3))**2 )
                    if (dij2 .lt. SurfCut2) then 
                      SoluteredoSA(k) = .true.
                      SolAtredoSA(k,l) = .true.
c                      goto 240
                    endif
 250              continue
 240            continue
              endif
c             Do Solvent--Proteins
              if (NProteins .gt. 0) then
                do 260 k=1,NProteins
                  nres = NResidues(k)
                  do 270 l=1,nres
c                    if (ResidueredoSA(k,l) .or. hasResidueMoved(k,l)) goto 270
                    resstart = BBZid(k,l,1)
                    resend = BBZid(k,l,4)
                    do 280 m=resstart,resend
                      dij2 = ( (xi-ZCoords(k,m,1))**2 + 
     &                  (yi-ZCoords(k,m,2))**2 +
     &                  (zi-ZCoords(k,m,3))**2 )
                      if (dij2 .lt. SurfCut2) then 
                        ResidueredoSA(k,l) = .true.
                        ProtAtredoSA(k,m) = .true.
c                        goto 270
                      endif
 280                 continue
 270               continue
 260             continue
              endif

 210         continue
          endif
 200     continue
      endif
c    
c     Checking what's scheduled for calc
c
c$$$      if (NSolutes .gt. 0) then 
c$$$        do 500 i=1,NSolutes
c$$$          if (IgnoreSolute(i)) goto 500
c$$$          write(printstring,*) 'Solute ',i,' doSA ',SoluteredoSA(i)
c$$$          call printLine(DEBUG,printstring)
c$$$          do 510 j=1,NSolAtoms(i)
c$$$            write (printstring,*) ' Solute ',i,' Atom ',j,' doSA ',SolAtredoSA(i,j)
c$$$            call printLine(DEBUG,printstring)
c$$$ 510      continue
c$$$ 500    continue
c$$$      endif
c$$$      
c$$$      if (NSolvents .gt. 0) then
c$$$        do 520 i=1,NSolvents
c$$$          write(printstring,*) 'Solvent ',i,' doSA ',SolventredoSA(i)
c$$$          call printLine(DEBUG,printstring)
c$$$          do 530 j=1,NSvnAtoms(i)
c$$$            write (printstring,*) ' Solvent ',i,' Atom ',j,' doSA ',SvnAtredoSA(i,j)
c$$$            call printLine(DEBUG,printstring)
c$$$ 530      continue
c$$$ 520    continue
c$$$      endif
c$$$
c$$$      if (NProteins .gt. 0) then
c$$$        do 540 i=1,NProteins
c$$$          do 550 j=1,NResidues(i)
c$$$            write(printstring,*) 'Prot ',i,' Res ',j,' doSA ',ResidueredoSA(i,j)
c$$$            call printLine(DEBUG,printstring)
c$$$            do 560 k=BBZid(i,j,1),BBZid(i,j,4)
c$$$              write (printstring,*) ' Protein ',i,' Residue ',j,
c$$$     &          ' Atom  ',k,' doSA',ProtAtredoSA(i,k)
c$$$              call printLine(DEBUG,printstring)
c$$$ 560        continue
c$$$ 550      continue
c$$$ 540    continue
c$$$      endif

      return
      end      
