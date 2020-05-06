      double precision function SolventSolventInterGBEnergy(isvn,jsvn,gbnrginterf,gbnrginterb)

      include 'dimensions.inc'
      include 'constants.inc'
      include 'solvents.inc'
      include 'printer.inc'
      include 'simparams.inc'
      include 'parameters.inc'
      include 'enums.inc'
      include 'flags.inc'
      include 'gb.inc'

      integer isvn,jsvn

      integer i,j,istart,iend,jstart,jend
c vactowat is -1 * ( 1 - 1 / 78.3 )
      double precision vactowat,distx,disty,distz,dist2,mindist,sclfac
      parameter (vactowat=-0.9872286079D+00)
      double precision gbnrginter,gbnrginterf,gbnrginterb
      double precision xi,yi,zi,borni,dij2,deno,bornij,temp
      
      logical imoved,jmoved
      integer istackdepth,jstackdepth
      double precision borniold,bornj,bornjold
      double precision xiold,yiold,ziold,dij2old
      double precision distxold,distyold,distzold,dist2old,mindistold
      double precision newgbnrg,oldgbnrg,newgbnrgf,newgbnrgb,oldgbnrgf,oldgbnrgb
      
      integer ipar,jpar
      logical ifepatom,jfepatom
      double precision bornif,bornib,bornioldf,bornioldb
      double precision qi,qi2,qif,qib,qif2,qib2
      double precision xif,yif,zif,xib,yib,zib
      double precision xioldf,yioldf,zioldf,xioldb,yioldb,zioldb
      double precision bornjf,bornjb,bornjoldf,bornjoldb
      double precision qj,qiqj,qiqjf,qiqjb
      double precision xj,yj,zj,xjf,yjf,zjf,xjb,yjb,zjb
      double precision xjold,yjold,zjold,xjoldf,yjoldf,zjoldf,xjoldb,yjoldb,zjoldb
      double precision dijf2,dijb2,dijf2old,dijb2old
      double precision bornijf,bornijb
      double precision bornijold,bornijfold,bornijbold

      gbnrginter = ZERO
      gbnrginterf = ZERO
      gbnrginterb = ZERO
      newgbnrg = ZERO
      oldgbnrg = ZERO
      newgbnrgf = ZERO
      oldgbnrgf = ZERO
      newgbnrgb = ZERO
      oldgbnrgb = ZERO
      mindist = 999999.0
      mindistold = 9999999.0


      imoved = .false.
      istackdepth = ZERO
      jmoved = .false.
      jstackdepth = ZERO

c
c Decide wether or not the interactions should be cut
c

      if (MovedSolvent(isvn) .gt. ZERO) then
        imoved = .true.
        istackdepth = MovedSolvent(isvn)
      endif
      if (MovedSolvent(jsvn) .gt. ZERO) then
        jmoved = .true.
        jstackdepth = MovedSolvent(jsvn)
      endif
      SolventSolventInterGBEnergy = gbnrginter
      if (imoved .and. jmoved) then
c       Old distance between the two residues
        distxold = dabs(gbSvnCentre(istackdepth,1)-gbSvnCentre(jstackdepth,1))
        distyold = dabs(gbSvnCentre(istackdepth,2)-gbSvnCentre(jstackdepth,2))
        distzold = dabs(gbSvnCentre(istackdepth,3)-gbSvnCentre(jstackdepth,3))
        dist2old = distxold**2 + distyold**2 + distzold**2
c       New distance between the two residues
        distx = dabs(SvnCentre(isvn,1)-SvnCentre(jsvn,1))
        disty = dabs(SvnCentre(isvn,2)-SvnCentre(jsvn,2))
        distz = dabs(SvnCentre(isvn,3)-SvnCentre(jsvn,3))
        dist2 = distx**2 + disty**2 + distz**2
        if (dist2old .gt. cut2 .and. dist2 .gt. cut2) return
      else if (imoved) then
c       Old distance between the two residues
        distxold = dabs(gbSvnCentre(istackdepth,1)-SvnCentre(jsvn,1))
        distyold = dabs(gbSvnCentre(istackdepth,2)-SvnCentre(jsvn,2))
        distzold = dabs(gbSvnCentre(istackdepth,3)-SvnCentre(jsvn,3))
        dist2old = distxold**2 + distyold**2 + distzold**2
c       New distance between the two residues
        distx = dabs(SvnCentre(isvn,1)-SvnCentre(jsvn,1))
        disty = dabs(SvnCentre(isvn,2)-SvnCentre(jsvn,2))
        distz = dabs(SvnCentre(isvn,3)-SvnCentre(jsvn,3))
        dist2 = distx**2 + disty**2 + distz**2
        if (dist2old .gt. cut2 .and. dist2 .gt. cut2) return        
      else if (jmoved) then
c       Old distance between the two residues
        distxold = dabs(SvnCentre(isvn,1)-gbSvnCentre(jstackdepth,1))
        distyold = dabs(SvnCentre(isvn,2)-gbSvnCentre(jstackdepth,2))
        distzold = dabs(SvnCentre(isvn,3)-gbSvnCentre(jstackdepth,3))
        dist2old = distxold**2 + distyold**2 + distzold**2
c       New distance between the two residues
        distx = dabs(SvnCentre(isvn,1)-SvnCentre(jsvn,1))
        disty = dabs(SvnCentre(isvn,2)-SvnCentre(jsvn,2))
        distz = dabs(SvnCentre(isvn,3)-SvnCentre(jsvn,3))
        dist2 = distx**2 + disty**2 + distz**2
        if (dist2old .gt. cut2 .and. dist2 .gt. cut2) return  
      else
        distx = dabs(SvnCentre(isvn,1)-SvnCentre(jsvn,1))
        disty = dabs(SvnCentre(isvn,2)-SvnCentre(jsvn,2))
        distz = dabs(SvnCentre(isvn,3)-SvnCentre(jsvn,3))
        dist2 = distx**2 + disty**2 + distz**2
        dist2old = dist2
        if (dist2 .gt. cut2) return  
      endif

      mindist = dist2
      mindistold = dist2old
c
c Select indexes
c

      istart = 1
      iend = NSvnAtoms(isvn)
      jstart = 1
      jend = NSvnAtoms(jsvn)

c      write(printstring,*) 'imoved ',imoved,' jmoved ',jmoved
c      call printLine(DEBUG,printstring)
      do 10 i=istart,iend
        ifepatom = .false.

        ipar = ZSvnPar(isvn,i)
        borni = ZSvnEffBorn(isvn,i)
        bornif = ZSvnEffBornF(isvn,i)
        bornib = ZSvnEffBornB(isvn,i)

        borniold = OldZSvnEffBorn(isvn,i)      
        bornioldf = OldZSvnEffBornF(isvn,i)      
        bornioldb = OldZSvnEffBornB(isvn,i)  

        if (ipar .lt. 0) then
          qi = VariableChg(-ipar,2)
          qi2 = qi ** 2
          if (DifferentLambdas) then
            ifepatom = .true.
            qif = VariableChg(-ipar,5)
            qib = VariableChg(-ipar,4)
            qif2 = qif ** 2
            qib2 = qib ** 2
          else
            qif = qi
            qib = qi
            qif2 = qi2
            qib2 = qi2
          endif
        else
          qi = ChgDB(CLJType(ipar))
          qi2 = qi ** 2
          qif = qi 
          qib = qi
          qif2 = qi2
          qib2 = qi2
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
        if (imoved) then
          xiold = gbSolventCoords(istackdepth,i,1)
          yiold = gbSolventCoords(istackdepth,i,2)
          ziold = gbSolventCoords(istackdepth,i,3)
          xioldf = xiold
          yioldf = yiold
          zioldf = ziold
          xioldb = xiold
          yioldb = yiold
          zioldb = ziold
        else
          xiold = xi
          yiold = yi
          ziold = zi
          xioldf = xif
          yioldf = yif
          zioldf = zif
          xioldb = xib
          yioldb = yib
          zioldb = zib
        endif

        do 20 j=jstart,jend
          jfepatom = .false.
c
c         See if this interaction needs to be computed
c

          if (.not. imoved .and. .not. jmoved .and.
     &        .not. SvnAtomdoBorn(isvn,i) .and. .not. SvnAtomdoBorn(jsvn,j)) then
            xj = ZSvnCoords(jsvn,j,1)
            yj = ZSvnCoords(jsvn,j,2)
            zj = ZSvnCoords(jsvn,j,3)

c         Still have to redo the distances for the cutoff
            dij2 = (xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2
            dij2old = dij2
c           Jump to the cutoff distance test
            goto 666
          endif

          jpar = ZSvnPar(jsvn,j)
          bornj = ZSvnEffBorn(jsvn,j)
          bornjf = ZSvnEffBornF(jsvn,j)
          bornjb = ZSvnEffBornB(jsvn,j)

          bornjold = OldZSvnEffBorn(jsvn,j)
          bornjoldf = OldZSvnEffBornF(jsvn,j)
          bornjoldb = OldZSvnEffBornB(jsvn,j)

          if (jpar .lt. 0) then
            qiqj = qi*VariableChg(-jpar,2)
            if (DifferentLambdas) then
              jfepatom = .true.
              qiqjf = qif*VariableChg(-jpar,5)
              qiqjb = qib*VariableChg(-jpar,4)
            else
              qiqjf = qiqj
              qiqjb = qiqj
            endif
          else
            qj = ChgDB(CLJType(jpar))
            qiqj = qi*qj
            qiqjf = qif*qj
            qiqjb = qib*qj
          endif

          xj = ZSvnCoords(jsvn,j,1)
          yj = ZSvnCoords(jsvn,j,2)
          zj = ZSvnCoords(jsvn,j,3)
          xjf = xj
          yjf = yj
          zjf = zj
          xjb = xj
          yjb = yj
          zjb = zj 
          if (jmoved) then
            xjold = gbSolventCoords(jstackdepth,j,1)
            yjold = gbSolventCoords(jstackdepth,j,2)
            zjold = gbSolventCoords(jstackdepth,j,3)
            xjoldf = xjold
            yjoldf = yjold
            zjoldf = zjold
            xjoldb = xjold
            yjoldb = yjold
            zjoldb = zjold
          else
            xjold = xj
            yjold = yj
            zjold = zj
            xjoldf = xjf
            yjoldf = yjf
            zjoldf = zjf
            xjoldb = xjb
            yjoldb = yjb
            zjoldb = zjb
          endif

C         Get the distances
          if (.not. ifepatom .and. .not. jfepatom) then
            dij2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
            dijf2 = dij2
            dijb2 = dij2    
            if (imoved .or. jmoved) then
              dij2old = (xiold-xjold)**2 + (yiold-yjold)**2 + (ziold-zjold)**2
              dijf2old = dij2old
              dijb2old = dij2old   
            else
              dij2old = dij2
              dijf2old = dij2
              dijb2old = dij2
            endif
          else
            dij2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
            dijf2 = (xif-xjf)**2 + (yif-yjf)**2 + (zif-zjf)**2
            dijb2 = (xib-xjb)**2 + (yib-yjb)**2 + (zib-zjb)**2            
            if (imoved .or. jmoved) then
              dij2old = (xiold-xjold)**2 + (yiold-yjold)**2 + (ziold-zjold)**2
              dijf2old = (xioldf-xjoldf)**2 + (yioldf-yjoldf)**2 + (zioldf-zjoldf)**2
              dijb2old = (xioldb-xjoldb)**2 + (yioldb-yjoldb)**2 + (zioldb-zjoldb)**2
            else
              dij2old = dij2
              dijf2old = dijf2
              dijb2old = dijb2
            endif
          endif
c         Get the product of the born radii
          bornij = borni*bornj
          bornijf = bornif*bornjf
          bornijb = bornib*bornjb          

          bornijold = borniold*bornjold
          bornijfold = bornioldf*bornjoldf
          bornijbold = bornioldb*bornjoldb
c         Compute the pairwise energy
c         Some optimisation is possible here

c         New energy -- reference state --
          deno = dsqrt(dij2 + bornij*dexp(-0.25d0*dij2/bornij))
          temp = qiqj / deno 
          newgbnrg = newgbnrg + temp
c         Old energy -- reference state --
          deno = dsqrt(dij2old + bornijold*dexp(-0.25d0*dij2old/bornijold))
          temp = qiqj / deno 
          oldgbnrg = oldgbnrg + temp
          if (DifferentLambdas) then
c           New Energy -- forwards state --
            deno = dsqrt(dijf2 + bornijf*dexp(-0.25d0*dijf2/bornijf))
            temp = qiqjf / deno 
            newgbnrgf = newgbnrgf + temp
c           Old Energy -- forwards state --
            deno = dsqrt(dijf2old + bornijfold*dexp(-0.25d0*dijf2old/bornijfold))
            temp = qiqjf / deno 
            oldgbnrgf = oldgbnrgf + temp
c           New Energy -- backwards state --
            deno = dsqrt(dijb2 + bornijb*dexp(-0.25d0*dijb2/bornijb))
            temp = qiqjb / deno 
            newgbnrgb = newgbnrgb + temp
c           Old Energy -- backwards state --
            deno = dsqrt(dijb2old + bornijbold*dexp(-0.25d0*dijb2old/bornijbold))
            temp = qiqjb / deno 
            oldgbnrgb = oldgbnrgb + temp
          else
            newgbnrgf = newgbnrg
            newgbnrgb = newgbnrg
            oldgbnrgf = oldgbnrg
            oldgbnrgb = oldgbnrg
          endif

 666      continue
          if (dij2 .lt. mindist) mindist = dij2
          if (dij2old .lt. mindistold) mindistold = dij2old
c          write(printstring,*) 'dij2',dij2
c          call printLine(DEBUG,printstring)
 20     continue
 10   continue

c
c Scale new/old differently according to the value of mindist /cut2/feather2
c
      if (mindist.lt.cut2) then
c       within cutoff distance
        if (mindist.gt.feather2) then
c         in the feather region
          sclfac = featherfac * (cut2-mindist)
          newgbnrg = newgbnrg * sclfac
          newgbnrgf = newgbnrgf * sclfac
          newgbnrgb = newgbnrgb * sclfac
        endif
      else
        newgbnrg = ZERO
        newgbnrgf = ZERO
        newgbnrgb = ZERO
      endif
      if (mindistold.lt.cut2) then
c       within cutoff distance
        if (mindistold.gt.feather2) then
c         in the feather region
          sclfac = featherfac * (cut2-mindistold)
          oldgbnrg = oldgbnrg * sclfac
          oldgbnrgf = oldgbnrgf * sclfac
          oldgbnrgb = oldgbnrgb * sclfac
        endif
      else
        oldgbnrg = ZERO
        oldgbnrgf = ZERO
        oldgbnrgb = ZERO
      endif
c     If initialisation, old should be ZERO
      if (.not. GotInitialTotalEnergy) then
        oldgbnrg = ZERO
        oldgbnrgf = ZERO
        oldgbnrgb = ZERO
      endif
      gbnrginter = (newgbnrg - oldgbnrg ) * vactowat * FAC4PIEPS0
      gbnrginterf = (newgbnrgf - oldgbnrgf ) * vactowat * FAC4PIEPS0
      gbnrginterb = (newgbnrgb - oldgbnrgb ) * vactowat * FAC4PIEPS0

c      write(printstring,*) 'NEWGBNRG ',newgbnrg*vactowat
c      call printLine(DEBUG,printstring)
c      write(printstring,*) 'OLDGBNRG ',oldgbnrg*vactowat
c     call printLine(DEBUG,printstring)

      SolventSolventInterGBEnergy = gbnrginter
c      write(printstring,*) 'isvn ',isvn,' jsvn ',jsvn,
c     &  ' EGBinter is ',SolventSolventInterGBEnergy
c      call printLine(DEBUG,printstring)

      return 
      end
