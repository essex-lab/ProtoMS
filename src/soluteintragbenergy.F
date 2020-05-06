      double precision function SoluteIntraGBEnergy(isol,gbnrgintraf,gbnrgintrab)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'solutes.inc'
      include 'printer.inc'
      include 'simparams.inc'
      include 'parameters.inc'
      include 'enums.inc'
      include 'flags.inc'
      include 'gb.inc'
      integer isol

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
      iend = NSolAtoms(isol)
      if (MovedSolute(isol) .gt. ZERO ) then
        moved = .true.
        stackdepth = MovedSolute(isol)
      endif


      do 10 i=istart,iend
        ifepatom = .false.

        ipar = ZSolPar(isol,i)
        borni = ZSolEffBorn(isol,i)
        bornif = ZSolEffBornF(isol,i)
        bornib = ZSolEffBornB(isol,i)

        borniold = OldZSolEffBorn(isol,i)      
        bornioldf = OldZSolEffBornF(isol,i)      
        bornioldb = OldZSolEffBornB(isol,i)      

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

        xi = ZSolCoords(isol,i,1)
        yi = ZSolCoords(isol,i,2)
        zi = ZSolCoords(isol,i,3)
        if (ifepatom .and. NSolVariable(isol) .gt. 0) then
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
        if (moved) then
          xiold = gbSoluteCoords(stackdepth,i,1)
          yiold = gbSoluteCoords(stackdepth,i,2)
          ziold = gbSoluteCoords(stackdepth,i,3)
          if (ifepatom .and. NSolVariable(isol) .gt. 0) then
            xioldf = gbSoluteCoordsF(stackdepth,i,1)
            yioldf = gbSoluteCoordsF(stackdepth,i,2)
            zioldf = gbSoluteCoordsF(stackdepth,i,3)  
            xioldb = gbSoluteCoordsB(stackdepth,i,1)
            yioldb = gbSoluteCoordsB(stackdepth,i,2)
            zioldb = gbSoluteCoordsB(stackdepth,i,3)  
          else
            xioldf = xiold
            yioldf = yiold
            zioldf = ziold
            xioldb = xiold
            yioldb = yiold
            zioldb = ziold              
          endif
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


        if (moved .or. SolAtomdoBorn(isol,i)) then

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
c        if (i .eq. 1) then
c         write(printstring,*) 'Sol ',isol,' at ',i,' qi ',qi,' qif ',qif,
c     &      ' qib ',qib,' borni ',borni,' bornif ',bornif,' bornib ',bornib,
c     &      ' borniold ',borniold,' bornioldf ',bornioldf,' bornioldb ',bornioldb
c          call printLine(DEBUG,printstring)
          
c        write(printstring,*) ' xi ',xi,' yi ',yi,' zi ',zi
c        call printLine(DEBUG,printstring)
c        write(printstring,*) ' xiold ',xiold,' yiold ',yiold,' ziold ',ziold
c        call printLine(DEBUG,printstring)
        
c          write(printstring,*) ' at ',i,
c     .      ' new gbnrg ',newgbnrg*vactowat*FAC4PIEPS0,
c     .      ' old gb ',oldgbnrg*vactowat*FAC4PIEPS0
c          call printLine(DEBUG,printstring)
c        endif
        do 20 j=i+1,iend
          jfepatom = .false.
C         Consider skipping interaction if molecule hasn't moved and the pair 
c         of born radii i and j hasn't changed too much

          if (.not. moved .and. .not. SolAtomdoBorn(isol,i) .and. 
     &      .not. SolAtomdoBorn(isol,j))  goto 20
          jpar = ZSolPar(isol,j)
          bornj = ZSolEffBorn(isol,j)
          bornjf = ZSolEffBornF(isol,j)
          bornjb = ZSolEffBornB(isol,j)

          bornjold = OldZSolEffBorn(isol,j)
          bornjoldf = OldZSolEffBornF(isol,j)
          bornjoldb = OldZSolEffBornB(isol,j)

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

          xj = ZSolCoords(isol,j,1)
          yj = ZSolCoords(isol,j,2)
          zj = ZSolCoords(isol,j,3)
          if (jfepatom .and. NSolVariable(isol) .gt. 0) then
            xjf = ZSolCoords(NSolForwards(isol),j,1)
            yjf = ZSolCoords(NSolForwards(isol),j,2)
            zjf = ZSolCoords(NSolForwards(isol),j,3)
            xjb = ZSolCoords(NSolBackwards(isol),j,1)
            yjb = ZSolCoords(NSolBackwards(isol),j,2)
            zjb = ZSolCoords(NSolBackwards(isol),j,3)
          else
            xjf = xj
            yjf = yj
            zjf = zj
            xjb = xj
            yjb = yj
            zjb = zj
          endif
          if (moved) then
            xjold = gbSoluteCoords(stackdepth,j,1)
            yjold = gbSoluteCoords(stackdepth,j,2)
            zjold = gbSoluteCoords(stackdepth,j,3)
            if (jfepatom .and. NSolVariable(isol) .gt. 0) then
              xjoldf = gbSoluteCoordsF(stackdepth,j,1)
              yjoldf = gbSoluteCoordsF(stackdepth,j,2)
              zjoldf = gbSoluteCoordsF(stackdepth,j,3)
              xjoldb = gbSoluteCoordsB(stackdepth,j,1)
              yjoldb = gbSoluteCoordsB(stackdepth,j,2)
              zjoldb = gbSoluteCoordsB(stackdepth,j,3)
            else
              xjoldf = xjold
              yjoldf = yjold
              zjoldf = zjold
              xjoldb = xjold
              yjoldb = yjold
              zjoldb = zjold              
            endif
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

c         write(printstring,*) ' xi ',xi,' xiold ',xiold,' xif ',xif,
c    .      ' xifold ',xioldf,' xib ',xib,' xibold ',xioldb
c         call printLine(DEBUG,printstring)
c         write(printstring,*) ' xj ',xj,' xjold ',xjold,' xjf ',xjf,
c    .      ' xjfold ',xjoldf,' xjb ',xjb,' xjbold ',xjoldb
c         call printLine(DEBUG,printstring)
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
c         write(printstring,*) 'i ',i,' borni ',borni,' bornif ',bornif,
c     .      ' bornib ',bornib
c         call printLine(DEBUG,printstring)
c         write(printstring,*) 'j ',j,' bornj ',bornj,' bornjf ',bornjf,
c     .      ' bornjb ',bornjb
c         call printLine(DEBUG,printstring)
c         write(printstring,*) ' i ',i,' j ',j,' dij2 ',dij2,' dij2old ',dij2old
c     .      ' dijf2 ',dijf2,' dijf2old ',dijf2old,
c     .      ' dijb2 ',dijb2,' dijb2old ',dijb2old
c         call printLine(DEBUG,printstring)
c         write(printstring,*) ' i ',i,' j ',j,' bornij ',bornij,' bornijold ',bornijold
c         call printLine(DEBUG,printstring)

c         write(printstring,*) 'qiqj ',qiqj,' qiqjf ',qiqjf,' qiqjb ',qiqjb
c         call printLine(DEBUG,printstring)
c         New energy -- reference state --
          deno = dsqrt(dij2 + bornij*dexp(-0.25d0*dij2/bornij))
          temp = qiqj / deno 
          newgbnrg = newgbnrg + temp
c          write(printstring,*) ' i ',i,' j ',j,' newgbnrg ',temp
c          call printLine(DEBUG,printstring)
c         Old energy -- reference state --
          deno = dsqrt(dij2old + bornijold*dexp(-0.25d0*dij2old/bornijold))
          temp = qiqj / deno 
          oldgbnrg = oldgbnrg + temp
c          write(printstring,*) ' i ',i,' j ',j,' oldgbnrg ',temp
c          call printLine(DEBUG,printstring)
          if (DifferentLambdas) then
c           New Energy -- forwards state --
            deno = dsqrt(dijf2 + bornijf*dexp(-0.25d0*dijf2/bornijf))
            temp = qiqjf / deno 
            newgbnrgf = newgbnrgf + temp
c           write(printstring,*) 'newgbnrgf ',temp
c           call printLine(DEBUG,printstring)
c           Old Energy -- forwards state --
            deno = dsqrt(dijf2old + bornijfold*dexp(-0.25d0*dijf2old/bornijfold))
            temp = qiqjf / deno 
            oldgbnrgf = oldgbnrgf + temp
c           write(printstring,*) 'oldgbnrg ',temp
c           call printLine(DEBUG,printstring)
c           New Energy -- backwards state --
            deno = dsqrt(dijb2 + bornijb*dexp(-0.25d0*dijb2/bornijb))
            temp = qiqjb / deno 
            newgbnrgb = newgbnrgb + temp
c           write(printstring,*) 'newgbnrgb ',temp
c           call printLine(DEBUG,printstring)
c           Old Energy -- backwards state --
            deno = dsqrt(dijb2old + bornijbold*dexp(-0.25d0*dijb2old/bornijbold))
            temp = qiqjb / deno 
            oldgbnrgb = oldgbnrgb + temp
c           write(printstring,*) 'oldgbnrg ',temp
c           call printLine(DEBUG,printstring)
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
c      write(printstring,*) 'NEWGBNRG ',newgbnrg*vactowat * FAC4PIEPS0
c      call printLine(DEBUG,printstring)
c      write(printstring,*) 'OLDGBNRG ',oldgbnrg*vactowat * FAC4PIEPS0
c      call printLine(DEBUG,printstring)
c     write(printstring,*) 'NEWGBNRGF ',newgbnrgf*vactowat * FAC4PIEPS0
c     call printLine(DEBUG,printstring)
c      write(printstring,*) 'OLDGBNRGF ',oldgbnrgf*vactowat * FAC4PIEPS0
c     call printLine(DEBUG,printstring)
c      write(printstring,*) 'NEWGBNRGB ',newgbnrgb*vactowat * FAC4PIEPS0
c      call printLine(DEBUG,printstring)
c      write(printstring,*) 'OLDGBNRGB ',oldgbnrgb*vactowat * FAC4PIEPS0
c      call printLine(DEBUG,printstring)

      SoluteIntraGBEnergy = gbnrgintra
c      write(printstring,*) 'solute ',isol,
c     &  ' EGBintra is ',SoluteIntraGBEnergy
c      call printLine(DEBUG,printstring)

      return
      end
