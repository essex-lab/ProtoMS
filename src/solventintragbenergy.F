      double precision function SolventIntraGBEnergy(isvn,gbnrgintraf,gbnrgintrab)

      include 'dimensions.inc'
      include 'constants.inc'
      include 'solvents.inc'
      include 'printer.inc'
      include 'simparams.inc'
      include 'parameters.inc'
      include 'enums.inc'
      include 'flags.inc'
      include 'gb.inc'

      integer isvn

      integer i,j,istart,iend
c vactowat is -1 * ( 1 - 1 / 78.3 )
      double precision vactowat
      parameter (vactowat=-0.9872286079D+00)
      double precision gbnrgintra,gbnrgintraf,gbnrgintrab
      double precision xi,yi,zi,borni,dij2,deno,bornij,temp

      double precision xiold,yiold,ziold

      double precision bornj,borniold,bornjold
      double precision newgbnrg,oldgbnrg,newgbnrgf,newgbnrgb,oldgbnrgf,oldgbnrgb
      logical moved
      integer stackdepth

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
      double precision dijf2,dijb2,dij2old,dijf2old,dijb2old
      double precision bornijf,bornijb
      double precision bornijold,bornijfold,bornijbold

      gbnrgintra = ZERO
      gbnrgintraf = ZERO
      gbnrgintrab = ZERO
      newgbnrg = ZERO
      oldgbnrg = ZERO
      newgbnrgf = ZERO
      oldgbnrgf = ZERO
      newgbnrgb = ZERO
      oldgbnrgb = ZERO
      moved = .false.
      stackdepth = ZERO


      istart = 1
      iend = NSvnAtoms(isvn)
      if (MovedSolvent(isvn) .gt. ZERO) then
        moved = .true.
        stackdepth = MovedSolvent(isvn)
      endif


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
        if (moved) then
          xiold = gbSolventCoords(stackdepth,i,1)
          yiold = gbSolventCoords(stackdepth,i,2)
          ziold = gbSolventCoords(stackdepth,i,3)
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


        if (moved .or. SvnAtomdoBorn(isvn,i)) then

C         Get the new self energy term --reference state --
          temp = (0.5*qi2/borni)
          newgbnrg = newgbnrg + temp
c         Get the old energy --reference state--
          temp = (0.5*qi2/borniold)
          oldgbnrg = oldgbnrg + temp
          if (DifferentLambdas) then
c           Get the new self energy term -- forwards state--
            temp = (0.5*qif2/bornif)
            newgbnrgf = newgbnrgf + temp
c           Get the old energy  -- forwards state--
            temp = (0.5*qif2/bornioldf)
            oldgbnrgf = oldgbnrgf + temp
c           Get the new self energy term -- backwards state--
            temp = (0.5*qib2/bornib)
            newgbnrgb = newgbnrgb + temp
c           Get the old energy  -- backwards state--
            temp = (0.5*qib2/bornioldb)
            oldgbnrgb = oldgbnrgb + temp
          else
            newgbnrgf = newgbnrg
            oldgbnrgf = oldgbnrg
            newgbnrgb = newgbnrg
            oldgbnrgb = oldgbnrg
          endif
        endif
c#if defined(RESIDUE)
c        write(printstring,*) 'Res ',ires,' at ',i,' qi ',qi,' qif ',qif,
c     &    ' qib ',qib,' borni ',borni,' bornif ',bornif,' bornib ',bornib,
c     &    ' borniold ',borniold,' bornioldf ',bornioldf,' bornioldb ',bornioldb
c        call printLine(DEBUG,printstring)
c#endif
        do 20 j=i+1,iend
          jfepatom = .false.
C         Consider skipping interaction if molecule hasn't moved and the pair 
c         of born radii i and j hasn't changed too much

          if (.not. moved .and. .not. SvnAtomdoBorn(isvn,i) .and. 
     &      .not.SvnAtomdoBorn(isvn,j))  goto 20
          jpar = ZSvnPar(isvn,j)
          bornj = ZSvnEffBorn(isvn,j)
          bornjf = ZSvnEffBornF(isvn,j)
          bornjb = ZSvnEffBornB(isvn,j)

          bornjold = OldZSvnEffBorn(isvn,j)
          bornjoldf = OldZSvnEffBornF(isvn,j)
          bornjoldb = OldZSvnEffBornB(isvn,j)

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

          xj = ZSvnCoords(isvn,j,1)
          yj = ZSvnCoords(isvn,j,2)
          zj = ZSvnCoords(isvn,j,3)
          xjf = xj
          yjf = yj
          zjf = zj
          xjb = xj
          yjb = yj
          zjb = zj 
          if (moved) then
            xjold = gbSolventCoords(stackdepth,j,1)
            yjold = gbSolventCoords(stackdepth,j,2)
            zjold = gbSolventCoords(stackdepth,j,3)
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
            if (moved) then
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
            if (moved) then
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
            oldgbnrgf = oldgbnrg
            newgbnrgb = newgbnrg
            oldgbnrgb = oldgbnrg            
          endif
 20     continue
 10   continue
c     If initialisation, old should be ZERO
      if (.not. GotInitialTotalEnergy) then
        oldgbnrg = ZERO
        oldgbnrgf = ZERO
        oldgbnrgb = ZERO
      endif
c needs to be multiplied by the dielectric constant difference
      gbnrgintra = (newgbnrg - oldgbnrg) * vactowat * FAC4PIEPS0
      gbnrgintraf = (newgbnrgf - oldgbnrgf) * vactowat * FAC4PIEPS0
      gbnrgintrab = (newgbnrgb - oldgbnrgb) * vactowat * FAC4PIEPS0
c      write(printstring,*) 'NEWGBNRG ',newgbnrg*vactowat
c      call printLine(DEBUG,printstring)
c      write(printstring,*) 'OLDGBNRG ',oldgbnrg*vactowat
c      call printLine(DEBUG,printstring)


      SolventIntraGBEnergy = gbnrgintra
c      write(printstring,*) 'svn ',isvn,
c     &  ' EGBintra is ',SolventIntraGBEnergy
c      call printLine(DEBUG,printstring)

      return
      end





