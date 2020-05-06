      subroutine SolventBornIntraPDA(savenrg,isvn)

c####################################################
c Computes the PDA Born intra Pair-Pair terms
c
c####################################################
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'flags.inc'
      include 'parameters.inc'
      include 'enums.inc'
      include 'gb.inc'

      integer isvn

      integer savenrg,nats,i,j,atstart,atend
      double precision borni,invborni,sborni,sborni2
      double precision bornj,invbornj,sbornj,sbornj2
      double precision xi,yi,zi,dij,invdij,cont
      double precision invlij,invlij2,invuij,invuij2
      double precision invlji,invlji2,invuji,invuji2      
      
      integer ipar,jpar
      double precision bornif,sbornif,bornib,sbornib
      logical ifepatom,jfepatom
      double precision invbornif,invbornib,sbornif2,sbornib2
      double precision xif,yif,zif,xib,yib,zib
      double precision xj,yj,zj,xjf,yjf,zjf,xjb,yjb,zjb
      double precision bornjf,sbornjf,bornjb,sbornjb
      double precision invbornjf,invbornjb,sbornjf2,sbornjb2
      double precision dijf,dijb
      double precision invdijf,invdijb
      double precision invlijf,invlijf2,invuijf,invuijf2
      double precision invljif,invljif2,invujif,invujif2      
      double precision invlijb,invlijb2,invuijb,invuijb2
      double precision invljib,invljib2,invujib,invujib2


      atstart = 1
      atend = NSvnAtoms(isvn)

      
c
c Do the intra object contributions
c
      nats = atend - atstart + 1
      if (nats .gt. 1) then
        do 10 i=atstart,atend-1
c         For each atom i...
          ifepatom = .false.

          ipar = ZSvnPar(isvn,i)

          if (ipar .lt. 0) then
            borni = VariableIntrBorn(-ipar,2)
            sborni = VariablePDAScale(-ipar,2)
            invborni = 1 / borni
            sborni2 = sborni ** 2
            if (DifferentLambdas) then
              ifepatom = .true.
              bornif = VariableIntrBorn(-ipar,5)
              sbornif = VariablePDAScale(-ipar,5)
              bornib = VariableIntrBorn(-ipar,4)
              sbornib = VariablePDAScale(-ipar,4)
              
              invbornif = 1 / bornif
              sbornif2 = sbornif ** 2
              invbornib = 1 / bornib
              sbornib2 = sbornib ** 2
            else
              bornif = borni
              sbornif = sborni
              invbornif = invborni
              sbornif2 = sborni2
              bornib = borni
              sbornib = sborni
              invbornib = invborni
              sbornib2 = sborni2
            endif
          else
            borni = IntrBornDB(CLJType(ipar))
            sborni = PDAScaleDB(CLJType(ipar))
            invborni = 1 / borni
            sborni2 = sborni ** 2
            bornif = borni
            sbornif = sborni
            invbornif = invborni
            sbornif2 = sborni2
            bornib = borni
            sbornib = sborni
            invbornib = invborni
            sbornib2 = sborni2            
          endif

          xi = ZSvnCoords(isvn,i,1)
          yi = ZSvnCoords(isvn,i,2)
          zi = ZSvnCoords(isvn,i,3)
          xif = xi
          yif = yi
          zif = zi
          xib = xi
          yib = yi
          zib = zi

          do 20 j=i+1,atend
            jfepatom = .false.

            jpar = ZSvnPar(isvn,j)

           if (jpar .lt. 0) then
              bornj = VariableIntrBorn(-jpar,2)
              sbornj = VariablePDAScale(-jpar,2)
              invbornj = 1 / bornj
              sbornj2 = sbornj**2
              if (DifferentLambdas) then
                jfepatom = .true.
                bornjf = VariableIntrBorn(-jpar,5)
                sbornjf = VariablePDAScale(-jpar,5)
                bornjb = VariableIntrBorn(-jpar,4)
                sbornjb = VariablePDAScale(-jpar,4)

                invbornjf = 1 / bornjf
                sbornjf2 = sbornjf ** 2
                invbornjb = 1 / bornjb
                sbornjb2 = sbornjb ** 2
              else
                bornjf = bornj
                sbornjf = sbornj
                invbornjf = invbornj
                sbornjf2 = sbornj2
                bornjb = bornj
                sbornjb = sbornj
                invbornjb = invbornj
                sbornjb2 = sbornj2
              endif
            else
              bornj = IntrBornDB(CLJType(jpar))
              sbornj = PDAScaleDB(CLJType(jpar))
              invbornj = 1 / bornj
              sbornj2 = sbornj ** 2
              bornjf = bornj
              sbornjf = sbornj
              invbornjf = invbornj
              sbornjf2 = sbornj2
              bornjb = bornj
              sbornjb = sbornj
              invbornjb = invbornj
              sbornjb2 = sbornj2
            endif

            xj = ZSvnCoords(isvn,j,1)
            yj = ZSvnCoords(isvn,j,2)
            zj = ZSvnCoords(isvn,j,3)
            xjf = xj
            yjf = yj
            zjf = zj
            xjb = xj
            yjb = yj
            zjb = zj      

            if (.not. ifepatom .and. .not. jfepatom) then
c             If the two atoms are non FEP then I just have to 
c             compute the contribution in the reference state and it will be the 
c             same in the forward/backward state
              dij = dsqrt ( (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2 )
              invdij = 1 / dij
              invuij = 1 / (dij + sbornj)
c             ...accumulate contributions from atom j
              if (invborni .gt. invuij) then
                invlij = 1 / dmax1(borni, dij - sbornj)
                invlij2 = invlij**2
                invuij2 = invuij**2
                cont = (invlij-invuij) + 
     &            0.25*dij*(invuij2-invlij2 ) + 
     &            0.5*invdij*dlog(invuij/invlij) +
     &            0.25*invdij*sbornj2*(invlij2-invuij2)

                if (savenrg .eq. NEW) then
                  NewZSvnContBorn(isvn,i) = NewZSvnContBorn(isvn,i) + cont
                  NewZSvnContBornF(isvn,i) = NewZSvnContBornF(isvn,i) + cont
                  NewZSvnContBornB(isvn,i) = NewZSvnContBornB(isvn,i) + cont
                else if (savenrg .eq. OLD) then
                  OldZSvnContBorn(isvn,i) = OldZSvnContBorn(isvn,i) + cont
                  OldZSvnContBornF(isvn,i) = OldZSvnContBornF(isvn,i) + cont
                  OldZSvnContBornB(isvn,i) = OldZSvnContBornB(isvn,i) + cont
                endif 

              endif
c             do contribution j--i as well
              invuji = 1 / (dij + sborni)
              if (invbornj .gt. invuji) then
                invlji = 1 / dmax1(bornj, dij - sborni)
                invlji2 = invlji**2
                invuji2 =invuji**2
                cont = (invlji-invuji) +
     &            0.25*dij*(invuji2-invlji2 ) + 
     &            0.5*invdij*dlog(invuji/invlji) +
     &            0.25*invdij*sborni2*(invlji2-invuji2)

                if (savenrg .eq. NEW) then
                  NewZSvnContBorn(isvn,j) = NewZSvnContBorn(isvn,j) + cont
                  NewZSvnContBornF(isvn,j) = NewZSvnContBornF(isvn,j) + cont
                  NewZSvnContBornB(isvn,j) = NewZSvnContBornB(isvn,j) + cont
                else if (savenrg .eq. OLD) then
                  OldZSvnContBorn(isvn,j) = OldZSvnContBorn(isvn,j) + cont
                  OldZSvnContBornF(isvn,j) = OldZSvnContBornF(isvn,j) + cont
                  OldZSvnContBornB(isvn,j) = OldZSvnContBornB(isvn,j) + cont
                endif

              endif
            else
c
c             One or both atoms are FEP and the reference,forward,backward
c             states must all be calculated 
c
              dij = dsqrt ( (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2 )
              invdij = 1 / dij
              invuij = 1 / (dij + sbornj)
c             ... accumulate contributions from atom j
              if (invborni .gt. invuij) then
                invlij = 1 / dmax1(borni, dij - sbornj)
                invlij2 = invlij**2
                invuij2 = invuij**2
                cont = (invlij-invuij) + 
     &            0.25*dij*(invuij2-invlij2 ) + 
     &            0.5*invdij*dlog(invuij/invlij) +
     &            0.25*invdij*sbornj2*(invlij2-invuij2)

                if (savenrg .eq. NEW) then
                  NewZSvnContBorn(isvn,i) = NewZSvnContBorn(isvn,i) + cont
                else if (savenrg .eq. OLD) then
                  OldZSvnContBorn(isvn,i) = OldZSvnContBorn(isvn,i) + cont
                endif   

              endif
c              ...do contribution j--i as well
              invuji = 1 / (dij + sborni)
              if (invbornj .gt. invuji) then
                invlji = 1 / dmax1(bornj, dij - sborni)
                invlji2 = invlji**2
                invuji2 =invuji**2
                cont = (invlji-invuji) +
     &            0.25*dij*(invuji2-invlji2 ) + 
     &            0.5*invdij*dlog(invuji/invlji) +
     &            0.25*invdij*sborni2*(invlji2-invuji2)

                if (savenrg .eq. NEW) then
                  NewZSvnContBorn(isvn,j) = NewZSvnContBorn(isvn,j) + cont
                else if (savenrg .eq. OLD) then
                  OldZSvnContBorn(isvn,j) = OldZSvnContBorn(isvn,j) + cont
                endif

              endif
c
c   Then do forwards calculation
c
              dijf = dsqrt ( (xif-xjf)**2 + (yif-yjf)**2 + (zif-zjf)**2 )
              invdijf = 1 / dijf
              invuijf = 1 / (dijf + sbornjf)
c             ... accumulate contributions from atom j
              if (invbornif .gt. invuijf) then
                invlijf = 1 / dmax1(bornif, dijf - sbornjf)
                invlijf2 = invlijf**2
                invuijf2 = invuijf**2
                cont = (invlijf-invuijf) + 
     &            0.25*dijf*(invuijf2-invlijf2 ) + 
     &            0.5*invdijf*dlog(invuijf/invlijf) +
     &            0.25*invdijf*sbornjf2*(invlijf2-invuijf2)      

                if (savenrg .eq. NEW) then
                  NewZSvnContBornF(isvn,i) = NewZSvnContBornF(isvn,i) + cont
                else if (savenrg .eq. OLD) then
                  OldZSvnContBornF(isvn,i) = OldZSvnContBornF(isvn,i) + cont
                endif   

              endif
c             do contribution from j--i
              invujif = 1 / (dijf + sbornif)
              if (invbornjf .gt. invujif) then
                invljif = 1 / dmax1(bornjf, dijf - sbornif)
                invljif2 = invljif**2
                invujif2 =invujif**2
                cont = (invljif-invujif) +
     &            0.25*dijf*(invujif2-invljif2 ) + 
     &            0.5*invdijf*dlog(invujif/invljif) +
     &            0.25*invdijf*sbornif2*(invljif2-invujif2)

                if (savenrg .eq. NEW) then
                  NewZSvnContBornF(isvn,j) = NewZSvnContBornF(isvn,j) + cont
                else if (savenrg .eq. OLD) then
                  OldZSvnContBornF(isvn,j) = OldZSvnContBornF(isvn,j) + cont
                endif   

              endif        
c
c  Then do backwards calculation
c
              dijb = dsqrt ( (xib-xjb)**2 + (yib-yjb)**2 + (zib-zjb)**2 )
              invdijb = 1 / dijb
              invuijb = 1 / (dijb + sbornjb)
c             ... accumulate contributions from atom j
              if (invbornib .gt. invuijb) then
                invlijb = 1 / dmax1(bornib, dijb - sbornjb)
                invlijb2 = invlijb**2
                invuijb2 = invuijb**2
                cont = (invlijb-invuijb) + 
     &            0.25*dijb*(invuijb2-invlijb2 ) + 
     &            0.5*invdijb*dlog(invuijb/invlijb) +
     &            0.25*invdijb*sbornjb2*(invlijb2-invuijb2)

                if (savenrg .eq. NEW) then
                  NewZSvnContBornB(isvn,i) = NewZSvnContBornB(isvn,i) + cont
                else if (savenrg .eq. OLD) then
                  OldZSvnContBornB(isvn,i) = OldZSvnContBornB(isvn,i) + cont
                endif   

              endif
c             do contribution j--i as well
              invujib = 1 / (dijb + sbornib)
              if (invbornjb .gt. invujib) then
                invljib = 1 / dmax1(bornjb, dijb - sbornib)
                invljib2 = invljib**2
                invujib2 =invujib**2
                cont = (invljib-invujib) +
     &            0.25*dijb*(invujib2-invljib2 ) + 
     &            0.5*invdijb*dlog(invujib/invljib) +
     &            0.25*invdijb*sbornib2*(invljib2-invujib2)

                if (savenrg .eq. NEW) then
                  NewZSvnContBornB(isvn,j) = NewZSvnContBornB(isvn,j) + cont
                else if (savenrg .eq. OLD) then
                  OldZSvnContBornB(isvn,j) = OldZSvnContBornB(isvn,j) + cont
                endif   

              endif
            endif
 20       continue
 10    continue    
      endif
      return 
      end
