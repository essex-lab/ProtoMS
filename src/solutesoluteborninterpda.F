      subroutine SoluteSoluteBornInterPDA(savenrg,isol,jsol)

      include 'dimensions.inc'
      include 'constants.inc'
      include 'simparams.inc'
      include 'solutes.inc'
      include 'flags.inc'
      include 'parameters.inc'
      include 'enums.inc'
      include 'gb.inc'

      integer isol,jsol

      integer savenrg,i,j,istart,iend,jstart,jend
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
      double precision scl,sclf,sclb
      double precision getLowScale,getHighScale

      if (isol .eq. jsol) return
      if (isol.le. 0 .or. isol.gt.MAXSOLUTES .or. jsol .le. 0 .or. jsol 
     &  .gt. MAXSOLUTES) return
      istart = 1
      iend = NSolAtoms(isol)
      jstart = 1
      jend = NSolAtoms(jsol)
      if (iend.le.0 .or. jend .le. 0) return  

      if (DualTopologySolute(isol).eq.jsol .or. DualTopologySolute(jsol).eq.isol) then
c       do not calculate the energy between dual topology pairs
c        write(printstring,*) 'sol ',isol,' and sol ',jsol,' NO NRG'
c       call printLine(DEBUG,printstring)
        return
      else if (DualTopologySolute(isol).eq.0 .and. DualTopologySolute(jsol).eq.0) then
c       there are no dual topology solutes
        scl = ONE
        sclb = ONE
        sclf = ONE
c       there are dual topology solutes, so get the scale factors
      else if (DualTopologySolute(isol).gt.0) then
c       solute 1 is a low dual topology solute
        scl = getLowScale(Lambda)
        sclb = getLowScale(LambdaB)
        sclf = getLowScale(LambdaF)
        if (DualTopologySolute(jsol).gt.0) then
c         solute 2 is also a low dual topology solute
          scl = scl*scl
          sclb = sclb*sclb
          sclf = sclf*sclf
        else if (DualTopologySolute(jsol).lt.0) then
c         solute 2 is a high dual topology solute
          scl = scl * getHighScale(Lambda)
          sclb = sclb * getHighScale(LambdaB)
          sclf = sclf * getHighScale(LambdaF)
        endif
      else if (DualTopologySolute(isol).eq.0) then
c       solute 1 is not a dual topology solute
        if (DualTopologySolute(jsol).gt.0) then
c         solute 2 is a low dual topology solute
          scl = getLowScale(Lambda)
          sclb = getLowScale(LambdaB)
          sclf = getLowScale(LambdaF)
        else
c         solute 2 is a high dual topology solute
          scl = getHighScale(Lambda)
          sclb = getHighScale(LambdaB)
          sclf = getHighScale(LambdaF)
        endif
      else
c       solute 1 is a high dual topology solute
        scl = getHighScale(Lambda)
        sclb = getHighScale(LambdaB)
        sclf = getHighScale(LambdaF)
        if (DualTopologySolute(jsol).gt.0) then
c         solute 2 is a low dual topology solute
          scl = scl * getLowScale(Lambda)
          sclb = sclb * getLowScale(LambdaB)
          sclf = sclf * getLowScale(LambdaF)
        else if (DualTopologySolute(jsol).lt.0) then
c         solute 2 is also a high dual topology solute
          scl = scl*scl
          sclb = sclb*sclb
          sclf = sclf*sclf
        endif
      endif

      do 10 i=istart,iend
c       For each atom i...
        ifepatom = .false.

        ipar = ZSolPar(isol,i)

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

        xi = ZSolCoords(isol,i,1)
        yi = ZSolCoords(isol,i,2)
        zi = ZSolCoords(isol,i,3)
        if (ifepatom .and. (NSolVariable(isol) .gt. 0)) then
          xif = ZSolCoords(NSolForwards(isol),i,1)
          yif = ZSolCoords(NSolForwards(isol),i,2)
          zif = ZSolCoords(NSolForwards(isol),i,3)
          xib = ZSolCoords(NSolBackwards(isol),i,1)
          yib = ZSolCoords(NSolBackwards(isol),i,2)
          zib = ZSolCoords(NSolBackwards(isol),i,3)
        else
          xif = xi
          yif = yi
          zif = zi
          xib = xi
          yib = yi
          zib = zi            
        endif

        do 20 j=jstart,jend
          jfepatom = .false.

          jpar = ZSolPar(jsol,j)

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

            xj = ZSolCoords(jsol,j,1)
            yj = ZSolCoords(jsol,j,2)
            zj = ZSolCoords(jsol,j,3)
            if (jfepatom .and. (NSolVariable(jsol) .gt. 0)) then
              xjf = ZSolCoords(NSolForwards(jsol),j,1)
              yjf = ZSolCoords(NSolForwards(jsol),j,2)
              zjf = ZSolCoords(NSolForwards(jsol),j,3)
              xjb = ZSolCoords(NSolBackwards(jsol),j,1)
              yjb = ZSolCoords(NSolBackwards(jsol),j,2)
              zjb = ZSolCoords(NSolBackwards(jsol),j,3)
            else
              xjf = xj
              yjf = yj
              zjf = zj
              xjb = xj
              yjb = yj
              zjb = zj
            endif

          if (.not. ifepatom .and. .not. jfepatom) then
c           If the two atoms are non FEP then I just have to 
c           compute the contribution in the reference state and it will be the 
c           same in the forward/backward state
            dij = dsqrt ( (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2 )
c           Consider cutting the interaction
            if (dij .gt. BornCut) goto 20
            invdij = 1 / dij
            invuij = 1 / (dij + sbornj)
c           ...accumulate contributions from atom j
            if (invborni .gt. invuij) then
              invlij = 1 / dmax1(borni, dij - sbornj)
              invlij2 = invlij**2
              invuij2 = invuij**2
              cont = (invlij-invuij) + 
     &          0.25*dij*(invuij2-invlij2 ) + 
     &          0.5*invdij*dlog(invuij/invlij) +
     &          0.25*invdij*sbornj2*(invlij2-invuij2)

              if (savenrg .eq. NEW) then
                NewZSolContBorn(isol,i) = NewZSolContBorn(isol,i) + cont * scl
                NewZSolContBornF(isol,i) = NewZSolContBornF(isol,i) + cont * sclf
                NewZSolContBornB(isol,i) = NewZSolContBornB(isol,i) + cont * sclb
              else if (savenrg .eq. OLD) then
                OldZSolContBorn(isol,i) = OldZSolContBorn(isol,i) + cont * scl
                OldZSolContBornF(isol,i) = OldZSolContBornF(isol,i) + cont * sclf
                OldZSolContBornB(isol,i) = OldZSolContBornB(isol,i) + cont * sclb
              endif                

            endif
c           do contribution j--i as well
            invuji = 1 / (dij + sborni)
            if (invbornj .gt. invuji) then
              invlji = 1 / dmax1(bornj, dij - sborni)
              invlji2 = invlji**2
              invuji2 =invuji**2
              cont = (invlji-invuji) +
     &          0.25*dij*(invuji2-invlji2 ) + 
     &          0.5*invdij*dlog(invuji/invlji) +
     &          0.25*invdij*sborni2*(invlji2-invuji2)

              if (savenrg .eq. NEW) then
                NewZSolContBorn(jsol,j) = NewZSolContBorn(jsol,j) + cont * scl
                NewZSolContBornF(jsol,j) = NewZSolContBornF(jsol,j) + cont * sclf
                NewZSolContBornB(jsol,j) = NewZSolContBornB(jsol,j) + cont * sclb
              else if (savenrg .eq. OLD) then
                OldZSolContBorn(jsol,j) = OldZSolContBorn(jsol,j) + cont * scl
                OldZSolContBornF(jsol,j) = OldZSolContBornF(jsol,j) + cont * sclf
                OldZSolContBornB(jsol,j) = OldZSolContBornB(jsol,j) + cont * sclb
              endif                

            endif
          else
c
c           One or both atoms are FEP and the reference,forward,backward
c           states must all be calculated 
c           
            dij = dsqrt ( (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2 )
c           Consider cutting the interaction
            if (dij .gt. BornCut) goto 20
            invdij = 1 / dij
            invuij = 1 / (dij + sbornj)
c           ... accumulate contributions from atom j
            if (invborni .gt. invuij) then
              invlij = 1 / dmax1(borni, dij - sbornj)
              invlij2 = invlij**2
              invuij2 = invuij**2
              cont = (invlij-invuij) + 
     &          0.25*dij*(invuij2-invlij2 ) + 
     &          0.5*invdij*dlog(invuij/invlij) +
     &          0.25*invdij*sbornj2*(invlij2-invuij2)

              if (savenrg .eq. NEW) then
                NewZSolContBorn(isol,i) = NewZSolContBorn(isol,i) + cont * scl
              else if (savenrg .eq. OLD) then
                OldZSolContBorn(isol,i) = OldZSolContBorn(isol,i) + cont * scl
              endif                

            endif              
c           ...do contribution j--i as well
            invuji = 1 / (dij + sborni)
            if (invbornj .gt. invuji) then
              invlji = 1 / dmax1(bornj, dij - sborni)
              invlji2 = invlji**2
              invuji2 =invuji**2
              cont = (invlji-invuji) +
     &          0.25*dij*(invuji2-invlji2 ) + 
     &          0.5*invdij*dlog(invuji/invlji) +
     &          0.25*invdij*sborni2*(invlji2-invuji2)            

              if (savenrg .eq. NEW) then
                NewZSolContBorn(jsol,j) = NewZSolContBorn(jsol,j) + cont * scl
              else if (savenrg .eq. OLD) then
                OldZSolContBorn(jsol,j) = OldZSolContBorn(jsol,j) + cont * scl
              endif                

            endif
 121        continue
c
c   Then do forwards calculation
c
            dijf = dsqrt ( (xif-xjf)**2 + (yif-yjf)**2 + (zif-zjf)**2 )
            if (dijf .gt. BornCut) goto 122
            invdijf = 1 / dijf
            invuijf = 1 / (dijf + sbornjf)
c           ... accumulate contributions from atom j
            if (invbornif .gt. invuijf) then
              invlijf = 1 / dmax1(bornif, dijf - sbornjf)
              invlijf2 = invlijf**2
              invuijf2 = invuijf**2
              cont = (invlijf-invuijf) + 
     &          0.25*dijf*(invuijf2-invlijf2 ) + 
     &          0.5*invdijf*dlog(invuijf/invlijf) +
     &          0.25*invdijf*sbornjf2*(invlijf2-invuijf2)                  

              if (savenrg .eq. NEW) then
                NewZSolContBornF(isol,i) = NewZSolContBornF(isol,i) + cont * sclf
              else if (savenrg .eq. OLD) then
                OldZSolContBornF(isol,i) = OldZSolContBornF(isol,i) + cont * sclf
              endif                

            endif                
c           do contribution from j--i
            invujif = 1 / (dijf + sbornif)
            if (invbornjf .gt. invujif) then
              invljif = 1 / dmax1(bornjf, dijf - sbornif)
              invljif2 = invljif**2
              invujif2 =invujif**2
              cont = (invljif-invujif) +
     &          0.25*dijf*(invujif2-invljif2 ) + 
     &          0.5*invdijf*dlog(invujif/invljif) +
     &          0.25*invdijf*sbornif2*(invljif2-invujif2)              

              if (savenrg .eq. NEW) then
                NewZSolContBornF(jsol,j) = NewZSolContBornF(jsol,j) + cont * sclf
              else if (savenrg .eq. OLD) then
                OldZSolContBornF(jsol,j) = OldZSolContBornF(jsol,j) + cont * sclf
              endif                

            endif   
 122        continue
c
c  Then do backwards calculation
c
            dijb = dsqrt ( (xib-xjb)**2 + (yib-yjb)**2 + (zib-zjb)**2 )
            if (dijb .gt. BornCut) goto 123
            invdijb = 1 / dijb
            invuijb = 1 / (dijb + sbornjb)
c           ... accumulate contributions from atom j
            if (invbornib .gt. invuijb) then
              invlijb = 1 / dmax1(bornib, dijb - sbornjb)
              invlijb2 = invlijb**2
              invuijb2 = invuijb**2
              cont = (invlijb-invuijb) + 
     &          0.25*dijb*(invuijb2-invlijb2 ) + 
     &          0.5*invdijb*dlog(invuijb/invlijb) +
     &          0.25*invdijb*sbornjb2*(invlijb2-invuijb2)            

              if (savenrg .eq. NEW) then
                NewZSolContBornB(isol,i) = NewZSolContBornB(isol,i) + cont * sclb
              else if (savenrg .eq. OLD) then
                OldZSolContBornB(isol,i) = OldZSolContBornB(isol,i) + cont * sclb
              endif                

            endif 
c           do contribution j--i as well
            invujib = 1 / (dijb + sbornib)
            if (invbornjb .gt. invujib) then
              invljib = 1 / dmax1(bornjb, dijb - sbornib)
              invljib2 = invljib**2
              invujib2 =invujib**2
              cont = (invljib-invujib) +
     &          0.25*dijb*(invujib2-invljib2 ) + 
     &          0.5*invdijb*dlog(invujib/invljib) +
     &          0.25*invdijb*sbornib2*(invljib2-invujib2)

              if (savenrg .eq. NEW) then
                NewZSolContBornB(jsol,j) = NewZSolContBornB(jsol,j) + cont * sclb
              else if (savenrg .eq. OLD) then
                OldZSolContBornB(jsol,j) = OldZSolContBornB(jsol,j) + cont * sclb
              endif                

            endif 
 123        continue
          endif
 20     continue
 10   continue    
      return 
      end 
