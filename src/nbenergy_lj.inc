c-------------INCLUDED FUNCTION AS AN INLINED BIT OF CODE-------------
c
c     This is only calculates the Coulombic part of the non-bonded energy
c
c     don't do anything if one of the parameters is null
      if (ipar1.ne.0.and.ipar2.ne.0) then
c***************************************************************
c       Work out the non-bonded parameters to use for this pair
c***************************************************************
        if (ipar1.lt.0 .or. ipar2.lt.0) then
c         interaction between variable clj parameters
c
c         funcsigsclf and funcsigsclb are used if one runs a simulation
c         with lorentz-berthelot combining rules and is reducing to zero
c         a sigma value (i.e creating a dummy atom in single topology)
c         This small modification tries to reduce the numerical instability 
c         associated with such perturbation.
c
          funcsigsclf = .false.
          funcsigsclb = .false.
          if (ipar1.lt.0) then           
            funcisig1 = VariableSig(-ipar1,2)
            funcieps1 = VariableEps(-ipar1,2)
            if (VariableSig(-ipar1,3).lt.TINY) funcsigsclf = .true.
            if (VariableSig(-ipar1,1).lt.TINY) funcsigsclb = .true.
          else
            funcisig1 = sigma(ipar1)
            funcieps1 = epslon(ipar1)
          endif
          if (ipar2.lt.0) then
            funcisig2 = VariableSig(-ipar2,2)
            funcieps2 = VariableEps(-ipar2,2)
            if (VariableSig(-ipar2,3).lt.TINY) funcsigsclf = .true.
            if (VariableSig(-ipar2,1).lt.TINY) funcsigsclb = .true.
          else
            funcisig2 = sigma(ipar2)
            funcieps2 = epslon(ipar2)
          endif
          if (lorentzberthelot) then
c           use lorentz-berthelot combining rules
c           (sigma is arithmetic mean, epsilon is geometric mean)
            funcsqrteps = FOUR * dsqrt(funcieps1*funcieps2)
            funcsig3 = (HALF*(funcisig1+funcisig2))**6
c           try to fix sigma scaling problem for disappearing atoms     
            if (funcsigsclf.and.funcsigsclb) then
              funcsig3 = ZERO
            else if (funcsigsclf) then
              funcsig3 = funcsig3 * (ONE-lambda)
            else if (funcsigsclb) then
              funcsig3 = funcsig3 * (lambda)
            endif
          else
c           use combining rules used in OPLS
c           (both sigma and epsilon are geometric means)
            funcsqrteps = FOUR * dsqrt(funcieps1*funcieps2)
            funcsig3 = (funcisig1*funcisig2)**3
          endif
c         *** Reference state non bonded parameters ***
          funcbij = funcsqrteps*funcsig3
          funcaij = funcbij*funcsig3
c         Calculate forwards and backwards free energies if there is a difference in 
c         lambda
          if (DifferentLambdas) then
            if (ipar1.lt.0) then               
               funcisig1b = VariableSig(-ipar1,4)
               funcisig1f = VariableSig(-ipar1,5)
               funcieps1b = VariableEps(-ipar1,4)
               funcieps1f = VariableEps(-ipar1,5)
            else              
               funcisig1f = funcisig1
               funcisig1b = funcisig1
               funcieps1f = funcieps1
               funcieps1b = funcieps1
            endif
            if (ipar2.lt.0) then
               funcisig2b = VariableSig(-ipar2,4)
               funcisig2f = VariableSig(-ipar2,5)
               funcieps2b = VariableEps(-ipar2,4)
               funcieps2f = VariableEps(-ipar2,5)
            else
               funcisig2f = funcisig2
               funcisig2b = funcisig2
               funcieps2f = funcieps2
               funcieps2b = funcieps2
            endif

c           forward and backward LJ energies
            if (lorentzberthelot) then
            
c             get the mixed sigma and epsilon parameters            
              funcsqrtepsf = FOUR * dsqrt(funcieps1f*funcieps2f)
              funcsig3f = ( HALF*(funcisig1f+funcisig2f) )**6
              funcsqrtepsb = FOUR * dsqrt(funcieps1b*funcieps2b)
              funcsig3b = ( HALF*(funcisig1b+funcisig2b) )**6

c             try to fix sigma scaling problem for disappering atoms 
              if (.not. SIRE_COMPATIBILITY_MODE ) then
                if (funcsigsclf.and.funcsigsclb) then
                  funcsig3f = ZERO
                  funcsig3b = ZERO
                else if (funcsigsclf) then
                  funcsig3f = funcsig3f * (ONE-lambdaf)
                  funcsig3b = funcsig3b * (ONE-lambdab)
                else if (funcsigsclb) then
                  funcsig3f = funcsig3f * (lambdaf)
                  funcsig3b = funcsig3b * (lambdab)
                endif
              endif
            else
              funcsqrtepsf = FOUR * dsqrt(funcieps1f*funcieps2f)
              funcsig3f = (funcisig1f*funcisig2f)**3
              funcsqrtepsb = FOUR * dsqrt(funcieps1b*funcieps2b)
              funcsig3b = (funcisig1b*funcisig2b)**3
            endif 
c           ***Forward/Backward state non bonded parameters ***         
            funcbijf = funcsqrtepsf*funcsig3f
            funcaijf = funcbijf*funcsig3f
            funcbijb = funcsqrtepsb*funcsig3b
            funcaijb = funcbijb*funcsig3b
          endif
        else
c         interaction between constant clj parameters
          funcbij = BIJ(CLJPair(ipar1,ipar2))
          funcaij = AIJ(CLJPair(ipar1,ipar2))
          funcbijf = funcbij
          funcbijb = funcbij
          funcaijf = funcaij
          funcaijb = funcaij
        endif
c**************************************************************
c       Work out the distances to use for this pair
c**************************************************************
c       if we are doing periodic boundaries then we need to find the nearest
c       image... - the higher level code (solutesoluteenergy.F etc.) has already
c       worked out whether we need to use periodic boundaries in the x/y/z dimensions
c       (signified by boundx/boundy/boundz)
        if (boundary.eq.PERIODICBOUNDARIES) then
          if (boundx) funcdx = funcdx - BoxDimension(1)
          if (boundy) funcdy = funcdy - BoxDimension(2)
          if (boundz) funcdz = funcdz - BoxDimension(3)
        endif

        funcdist2 = funcdx**2+funcdy**2+funcdz**2
        funcmindist = funcdist2
c       calculate the reciprocal of the distance here... 
c        funcrealdist = funcdist2
c       The softcore implemented here is the 
c       Separation shifted scaling soft core of Zacharias and McCammon
c       (reference here)
        if (softcore) then 
c         Check if we have alread worked out sigmaij6...
          if (ipar1 .lt. 0 .or. ipar2 .lt.0) then
c           Calculate the shift parameter
c           Taken as the sigma value for this pair pair interaction
           if (lorentzberthelot) then 
             shift = HALF*(funcisig1+funcisig2)
           else
             shift = funcisig1*funcisig2
           endif
          else
c           ...Otherwise we have to get it by lookup
            sigat1 = sigma(ipar1)
            sigat2 = sigma(ipar2)
            if (lorentzberthelot) then
              shift = HALF*(sigat1+sigat2)
            else
              shift = ( sigat1*sigat2 )
            endif
          endif
          if (SoftType.eq.SOFT66) then
            softdist = (1-scl)*DeltaCore*shift**6.0d0
            invfuncdist6 = ONE / (funcdist2**3 + softdist)
          else if (SoftType.eq.SOFTAMB) then
            softdist = (1-scl)*DeltaCore*shift**6.0d0
            invfuncdist6 = ONE / (funcdist2**3 + softdist)
          else
            softdist = (1-scl)*DeltaCore*shift
            invfuncdist6 = ONE / (funcdist2 + softdist)**3
          endif
          invfuncdist12 = invfuncdist6**2
        else
          invfuncdist6 = ONE / funcdist2**3
          invfuncdist12 = invfuncdist6**2
        endif
c       calculate the distances for variable geometry (if necessary)
        if (funcdogeom .and. DifferentLambdas) then
          if (boundx) then
            funcdxf = funcdxf - BoxDimension(1)
            funcdxb = funcdxb - BoxDimension(1)
          endif
          if (boundy) then
            funcdyf = funcdyf - BoxDimension(2)
            funcdyb = funcdyb - BoxDimension(2)
          endif
          if (boundz) then
            funcdzf = funcdzf - BoxDimension(3)
            funcdzb = funcdzb - BoxDimension(3)
          endif
          funcdistf2 = funcdxf**2+funcdyf**2+funcdzf**2
          funcdistb2 = funcdxb**2+funcdyb**2+funcdzb**2
          if (softcore) then            
            if (ipar1 .lt. 0 .or. ipar2 .lt.0) then
              if (lorentzberthelot) then 
                shiftf = HALF*(funcisig1f+funcisig2f)
                shiftb = HALF*(funcisig1b+funcisig2b)
              else
                shiftf = funcisig1f*funcisig2f
                shiftb = funcisig1b*funcisig2b
              endif
            else
c             ...Otherwise we have to get it by lookup
              if (lorentzberthelot) then
                shiftf =  HALF*(sigat1+sigat2 )
                shiftb =  shiftf
              else
                shiftf = ( sigat1*sigat2 )
                shiftb = shiftf
              endif
            endif
            if (SoftType.eq.SOFT66) then
              softdistf = (1-sclf)*DeltaCore*shiftf**6.0d0
              invfuncdist6f = ONE / (funcdistf2**3 + softdistf)
              softdistb = (1-sclb)*DeltaCore*shiftb**6.0d0
              invfuncdist6b = ONE / (funcdistb2**3 + softdistb)
            else if (SoftType.eq.SOFTAMB) then
              softdistf = (1-sclf)*DeltaCore*shiftf**6.0d0
              invfuncdist6f = ONE / (funcdistf2**3 + softdistf)
              softdistb = (1-sclb)*DeltaCore*shiftb**6.0d0
              invfuncdist6b = ONE / (funcdistb2**3 + softdistb)
            else
              softdistf = (1-sclf)*DeltaCore*shiftf
              invfuncdist6f = ONE / (funcdistf2 + softdistf)**3
              softdistb = (1-sclb)*DeltaCore*shiftb
              invfuncdist6b = ONE / (funcdistb2 + softdistb)**3
            endif
            invfuncdist12f = invfuncdist6f**2        
            invfuncdist12b = invfuncdist6b**2               
          else
            invfuncdist6f = ONE / funcdistf2**3
            invfuncdist12f = invfuncdist6f**2

            invfuncdist6b = ONE / funcdistb2**3
            invfuncdist12b = invfuncdist6b**2
          endif
        else
          if (softcore) then            
            if (ipar1 .lt. 0 .or. ipar2 .lt.0) then
              if (lorentzberthelot) then 
                shiftf = HALF*(funcisig1f+funcisig2f)
                shiftb = HALF*(funcisig1b+funcisig2b)
              else
                shiftf = funcisig1f*funcisig2f
                shiftb = funcisig1b*funcisig2b
              endif
            else
c             ...Otherwise we have to get it by lookup
              if (lorentzberthelot) then
                shiftf =  HALF*(sigat1+sigat2 )
                shiftb =  shiftf
              else
                shiftf = ( sigat1*sigat2 )
                shiftb = shiftf
              endif
            endif
            if (SoftType.eq.SOFT66) then
              softdistf = (1-sclf)*DeltaCore*shiftf**6.0d0
              invfuncdist6f = ONE / (funcdist2**3 + softdistf)
              softdistb = (1-sclb)*DeltaCore*shiftb**6.0d0
              invfuncdist6b = ONE / (funcdist2**3+ softdistb)
            else if (SoftType.eq.SOFTAMB) then
              softdistf = (1-sclf)*DeltaCore*shiftf**6.0d0
              invfuncdist6f = ONE / (funcdist2**3 + softdistf)
              softdistb = (1-sclb)*DeltaCore*shiftb**6.0d0
              invfuncdist6b = ONE / (funcdist2**3 + softdistb)
            else
              softdistf = (1-sclf)*DeltaCore*shiftf
              invfuncdist6f = ONE / (funcdist2 + softdistf)**3
              softdistb = (1-sclb)*DeltaCore*shiftb
              invfuncdist6b = ONE / (funcdist2 + softdistb)**3
            endif
            invfuncdist12f = invfuncdist6f**2        
            invfuncdist12b = invfuncdist6b**2        
          else
            invfuncdist6f = invfuncdist6
            invfuncdist6b = invfuncdist6
            invfuncdist12f = invfuncdist12
            invfuncdist12b = invfuncdist12
          endif
        endif
c***********************************************************
c       Now calculate the energies
c***********************************************************
c        call printLine(DEBUG,printstring)
        funccnrg = ZERO
        funcljnrg = funcaij*invfuncdist12 - funcbij*invfuncdist6
c        write(6,*)"f(r)= ",ONE/invfuncdist6**(ONE/6.0d0), "g(r)= ",ONE/invfuncdist
        if (DifferentLambdas) then
c         *** Forwards/Backwards state coulombic and LJ energies *** 
          funccnrgf = ZERO
          funcljnrgf = funcaijf*invfuncdist12f - funcbijf*invfuncdist6f
          funccnrgb = ZERO
          funcljnrgb = funcaijb*invfuncdist12b - funcbijb*invfuncdist6b
c          write(6,*)"f(r)= ",ONE/invfuncdist6f**(ONE/6.0d0), "g(r)= ",ONE/invfuncdistf
c          write(6,*)"f(r)= ",ONE/invfuncdist6b**(ONE/6.0d0), "g(r)= ",ONE/invfuncdistb
        else
          funccnrgf = ZERO
          funcljnrgf = funcljnrg
          funccnrgb = ZERO
          funcljnrgb = funcljnrg 
        endif
      endif
c      funcmindist = funcdist2
      funcnbnrg = funccnrg + funcljnrg
      funcnbnrgf = funccnrgf + funcljnrgf
      funcnbnrgb = funccnrgb + funcljnrgb

c     you can now safely read the values from 
c     funcnbnrg   =  non-bonded energy
c     funccnrg = coulomb component
c     funcljnrg = lj component
c     funcdist = distance between atoms

c     if there are different lambdas, then you will also find...
c     funcnbnrgf = forwards non-bonded energy
c     funcnbnrgb = backwards non-bonded energy
c     funccnrgf/funccnrgb = forwards/backwards coulomb energies
c     funcljnrgf/funcljnrgb = for/back LJ energies

c--------------END OF INCLUDED FUNCTION----------------------------------



    

