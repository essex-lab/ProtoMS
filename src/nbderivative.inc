c     don't do anything if one of the parameters is null
      if (ipar1.ne.0.and.ipar2.ne.0) then
c***************************************************************
c       Work out the non-bonded parameters to use for this pair
c***************************************************************
        if (ipar1.lt.0 .or. ipar2.lt.0) then
c         interaction between variable clj parameters
c         this is not supported, so set everything to zero
          qiqj = ZERO
          funcbij = ZERO
          funcaij = ZERO
          qiqjf = ZERO
          qiqjb = ZERO
          funcbijf = ZERO
          funcbijb = ZERO
          funcaijf = ZERO
          funcaijb = ZERO 
        else
c         interaction between constant clj parameters
          qiqj = ChgIJ(CLJPair(ipar1,ipar2))
          funcbij = BIJ(CLJPair(ipar1,ipar2))
          funcaij = AIJ(CLJPair(ipar1,ipar2))
          qiqjf = qiqj
          qiqjb = qiqj
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
        funcdist6 = funcdist2*funcdist2*funcdist2
        funcmindist = funcdist2
        if (softcore) then 
c         Calculate Coulomb part
          
          if (SoftType.eq.SOFT66) then
            softdist = funcdist6 + (1-scl)*DeltaCoreCoul
            invfuncdist = ONE / (softdist)**SoftPowerInv
            funccnrg = - CoulSoftPower*scl**(CoulSoftPower-1.0)*qiqj * invfuncdist 
     .                  - scl**CoulSoftPower*qiqj*DeltaCoreCoul * invfuncdist / (6.0d0*softdist)
          else if (SoftType.eq.SOFTAMB) then
            softdist = funcdist2 + (1-scl)*DeltaCoreCoul 
            invfuncdist = ONE / dsqrt( softdist )
            funccnrg = - CoulSoftPower*scl**(CoulSoftPower-1.0)*qiqj * invfuncdist -
     .                  scl**CoulSoftPower*qiqj*DeltaCoreCoul * invfuncdist / (2.0d0*softdist)
          else
            softdist = funcdist2 + (1-scl) 
            invfuncdist = ONE / dsqrt( softdist )
            funccnrg = - CoulSoftPower*scl**(CoulSoftPower-1.0)*qiqj * invfuncdist -
     .                  scl**CoulSoftPower*qiqj * invfuncdist / (2.0d0*softdist)
          endif
c         Check if we have alread worked out sigmaij6...
          if (ipar1 .lt. 0 .or. ipar2 .lt.0) then
c           this is not supported, so set it to zero
            shift = ZERO
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
c         Calculate Lennard-Jone part
          if ((SoftType.eq.SOFT66).or.(SoftType.eq.SOFTAMB)) then
            shift = DeltaCore*shift**6.0d0
            softdist = (1-scl)*shift+funcdist6
            softdist2 = softdist*softdist
            funcljnrg = funcbij / softdist - funcaij / softdist2 +
     .                  scl*funcbij*shift / softdist2 - scl*2.0d0*funcaij*shift / (softdist*softdist2)
          else
            shift = DeltaCore*shift
            softdist = (1-scl)*shift+funcdist2
            softdist2 = softdist*softdist*softdist
            softdist3 = softdist2*softdist2
            funcljnrg = funcbij / softdist2 - funcaij / softdist3 +
     .                  scl*3.0d0*funcbij*shift / (softdist*softdist2) - scl*6.0d0*funcaij*shift / (softdist*softdist3)
          endif
          invfuncdist12 = invfuncdist6**2
        else
          invfuncdist = ONE / dsqrt(funcdist2)
          invfuncdist6 = invfuncdist**6
          invfuncdist12 = invfuncdist6**2
          funccnrg = - CoulSoftPower*scl**(CoulSoftPower-1.0)*qiqj * invfuncdist
          funcljnrg = funcbij * invfuncdist6 - funcaij * invfuncdist12
        endif
c       Calculating derivatives at different lambdas make no sense
        funccnrgf = funccnrg
        funcljnrgf = funcljnrg
        funccnrgb = funccnrg
        funcljnrgb = funcljnrg 
      endif
      if (.not.doLowDerivative) then
        funccnrg = - funccnrg
        funcljnrg = - funcljnrg
      endif
      funcnbnrg = funccnrg + funcljnrg
      funcnbnrgf = funccnrgf + funcljnrgf
      funcnbnrgb = funccnrgb + funcljnrgb

c--------------END OF INCLUDED FUNCTION----------------------------------



    

