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

          if (ipar1.lt.0) then
            if (pertchg) then
              funcichg1 = VariableChg(-ipar1,2)
            else
              funcichg1 = VariableChg(-ipar1,1)
            endif            
          else
            funcichg1 = charge(ipar1)
          endif
          if (ipar2.lt.0) then
            if (pertchg) then
              funcichg2 = VariableChg(-ipar2,2)
            else
              funcichg2 = VariableChg(-ipar2,1)
            endif
          else
            funcichg2 = charge(ipar2)
          endif
c         *** Reference state non bonded parameters ***
          qiqj = funcichg1*funcichg2*FAC4PIEPS0
c         Calculate forwards and backwards free energies if there is a difference in 
c         lambda
          if (DifferentLambdas) then
            if (ipar1.lt.0) then
               funcichg1b = VariableChg(-ipar1,4)
               funcichg1f = VariableChg(-ipar1,5)
            else
               funcichg1b = funcichg1
               funcichg1f = funcichg1
            endif
            if (ipar2.lt.0) then
               funcichg2b = VariableChg(-ipar2,4)
               funcichg2f = VariableChg(-ipar2,5)
            else
               funcichg2b = funcichg2
               funcichg2f = funcichg2
            endif 
c           This turns off perturbation of charges
            if (.not.pertchg) then
              funcichg1b = funcichg1
              funcichg1f = funcichg1  
              funcichg2b = funcichg2
              funcichg2f = funcichg2  
            endif

c           ***Forward/Backward state non bonded parameters ***
            qiqjf = funcichg1f*funcichg2f*FAC4PIEPS0           
            qiqjb = funcichg1b*funcichg2b*FAC4PIEPS0
          endif
        else
c         interaction between constant clj parameters
          qiqj = ChgIJ(CLJPair(ipar1,ipar2))        
          qiqjf = qiqj
          qiqjb = qiqj
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
        if (softcore) then 
c         Soften the distance for the coulombic term
          if (SoftType.eq.SOFT66) then
            funcdist6 = funcdist2*funcdist2*funcdist2
            invfuncdist = ONE / (funcdist6 + (1-scl)*DeltaCoreCoul)**SoftPowerInv
          else if (SoftType.eq.SOFTAMB) then
            invfuncdist = ONE / dsqrt( funcdist2 + (1-scl)*DeltaCoreCoul )
          else
            invfuncdist = ONE / dsqrt( funcdist2 + (1-scl) )
          endif
        else
          invfuncdist = ONE / dsqrt(funcdist2)
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
            if (SoftType.eq.SOFT66) then
              funcdistf6 = funcdistf2*funcdistf2*funcdistf2
              invfuncdistf = ONE / (funcdistf6 + (1-sclf)*DeltaCoreCoul)**SoftPowerInv
              funcdistb6 = funcdistb2*funcdistb2*funcdistb2
              invfuncdistb = ONE / (funcdistb6 + (1-sclb)*DeltaCoreCoul)**SoftPowerInv
            else if (SoftType.eq.SOFTAMB) then
              invfuncdistf = ONE / dsqrt( funcdistf2 + (1-sclf)*DeltaCoreCoul )
              invfuncdistb = ONE / dsqrt( funcdistb2 + (1-sclb)*DeltaCoreCoul )
            else
              invfuncdistf = ONE / dsqrt( funcdistf2 + (1-sclf) )
              invfuncdistb = ONE / dsqrt( funcdistb2 + (1-sclb) )
            endif                     
          else
            invfuncdistf = ONE / dsqrt(funcdistf2)
            invfuncdistb = ONE / dsqrt(funcdistb2)
          endif
        else
          if (softcore) then
            if (SoftType.eq.SOFT66) then              
              invfuncdistf = ONE / (funcdist6 + (1-sclf)*DeltaCoreCoul)**SoftPowerInv
              invfuncdistb = ONE / (funcdist6 + (1-sclb)*DeltaCoreCoul)**SoftPowerInv
            else if (SoftType.eq.SOFTAMB) then
              invfuncdistf = ONE / dsqrt( funcdist2 + (1-sclf)*DeltaCoreCoul )
              invfuncdistb = ONE / dsqrt( funcdist2 + (1-sclb)*DeltaCoreCoul )
            else
              invfuncdistf = ONE / dsqrt( funcdist2 + (1-sclf) )
              invfuncdistb = ONE / dsqrt( funcdist2 + (1-sclb) )
            endif            
          else
            invfuncdistf = invfuncdist
            invfuncdistb = invfuncdist
          endif
        endif
c***********************************************************
c       Now calculate the energies
c***********************************************************
        if (pmeblocks.gt.1) then
          funccnrg = qiqj * invfuncdist * erfc(pme_alpha/invfuncdist)
        else
          funccnrg = qiqj * invfuncdist
        endif  
        funcljnrg = ZERO
        if (DifferentLambdas) then
c         *** Forwards/Backwards state coulombic and LJ energies *** 
          if (pmeblocks.gt.1) then
            funccnrgf = qiqjf * invfuncdistf * erfc(pme_alpha/invfuncdistf)
            funccnrgb = qiqjb * invfuncdistb * erfc(pme_alpha/invfuncdistb)
          else
            funccnrgf = qiqjf * invfuncdistf
            funccnrgb = qiqjb * invfuncdistb
          endif
          funcljnrgf = ZERO
          funcljnrgb = ZERO
        else
          funccnrgf = funccnrg
          funcljnrgf = ZERO
          funccnrgb = funccnrg
          funcljnrgb = ZERO
        endif
        if (ddd) then
          funccnrg = funccnrg * (0.25d0*invfuncdist)
          if (DifferentLambdas) then 
            funccnrgf = funccnrgf * (0.25d0*invfuncdistf)
            funccnrgb = funccnrgb * (0.25d0*invfuncdistb)
          else
            funccnrgf = funccnrg
            funccnrgb = funccnrg
          endif
c *********************************************************************
c ***************** Sigmoidal DDD function ****************************
c$$$          fac = -20.929 + 99.329 / ( 1 + 3.4781 * 
c$$$     .      dexp(-0.001787 * 99.329 * ( ONE / invfuncdist )))
c$$$c          write(printstring,*) 'Fac ',fac,' funccnrg was ',funccnrg, 
c$$$c     .      ' dist is ',ONE / invfuncdist
c$$$c          call printLine(DEBUG,printstring)
c$$$          funccnrg = funccnrg / fac
c$$$c          write(printstring,*) ' funccnrg is now ',funccnrg
c$$$c          call printLine(DEBUG,printstring)
c$$$          if (DifferentLambdas) then
c$$$            facf = -20.929 + 99.329 / ( 1 + 3.4781 * 
c$$$     .      dexp(-0.001787 * 99.329 * ( ONE / invfuncdistf )))
c$$$            funccnrgf = funccnrgf / facf
c$$$            facb = -20.929 + 99.329 / ( 1 + 3.4781 * 
c$$$     .      dexp(-0.001787 * 99.329 * ( ONE / invfuncdistb )))
c$$$            funccnrgb = funccnrgb / facb
c$$$          else
c$$$            funccnrgf = funccnrg
c$$$            funccnrgb = funccnrg
c$$$          endif
c *********************************************************************
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



    

