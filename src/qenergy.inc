c-------------INCLUDED FUNCTION AS AN INLINED BIT OF CODE-------------
c-------------This is because the function call was too expensive!----
c
c This is a rewrite of nbenergy2.inc to only calculate Coulomb interaction
c the charge initialisation is done by qinitial.inc
c     don't do anything if one of the parameters is null
      if (ipar1.ne.0.and.ipar2.ne.0) then
c**************************************************************
c       Work out the distances to use for this pair
c**************************************************************
        if (boundary.eq.PERIODICBOUNDARIES) then
          if (boundx) funcdx = funcdx - BoxDimension(1)
          if (boundy) funcdy = funcdy - BoxDimension(2)
          if (boundz) funcdz = funcdz - BoxDimension(3)
        endif

        funcdist2 = funcdx**2+funcdy**2+funcdz**2
c       calculate the reciprocal of the distance here... 
c        funcrealdist = funcdist2
c       The softcore implemented here is the 
c       Separation shifted scaling soft core of Zacharias and McCammon
c       (reference here)
        if (softcore) then 
c         Soften the distance for the coulombic term
          invfuncdist = ONE / dsqrt( funcdist2 + (1-scl) )
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
            invfuncdistf = ONE / dsqrt(funcdistf2 + sclf )
            invfuncdistb = ONE / dsqrt(funcdistb2 + sclb)                    
          else
            invfuncdistf = ONE / dsqrt(funcdistf2)
            invfuncdistb = ONE / dsqrt(funcdistb2)
          endif
        else
          if (softcore) then
            invfuncdistf = ONE / dsqrt(funcdist2 + (1-sclf) )
            invfuncdistb = ONE / dsqrt(funcdist2 + (1-sclb) )          
          else
            invfuncdistf = invfuncdist
            invfuncdistb = invfuncdist
          endif
        endif
c***********************************************************
c       Now calculate the energies
c***********************************************************
        funccnrg = qiqj * invfuncdist * erf(pme_alpha/invfuncdist)      
        if (DifferentLambdas) then
c         *** Forwards/Backwards state coulombic and LJ energies *** 
          funccnrgf = qiqjf * invfuncdistf * erf(pme_alpha/invfuncdistf)
          funccnrgb = qiqjb * invfuncdistb * erf(pme_alpha/invfuncdistb)
        else
          funccnrgf = funccnrg
          funccnrgb = funccnrg
        endif
      
      endif

      funcmindist = funcdist2
      funcnbnrg = funccnrg
      funcnbnrgf = funccnrgf
      funcnbnrgb = funccnrgb

c--------------END OF INCLUDED FUNCTION----------------------------------



    

