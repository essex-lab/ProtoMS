c-------------INCLUDED FUNCTION AS AN INLINED BIT OF CODE-------------
c-------------This is because the function call was too expensive!----
c
c This is a rewrite of nbenergy2.inc so that soft core calculations 
c are possible
c
c-- This function calculates the non-bonded energy between atoms
c-- whose coordinates are placed in at1(3) and at2(3), and whose
c-- parameters are in ipar1, and ipar2. The returned energy is called
c-- nbEnergy, the coulomb component is in funccnrg, the LJ component
c-- is in funcljnrg, and the distance is is funcmindist.
c
C     To use this function put the type of system into functyp
c     (PROTEIN, SOLUTE or SOLVENT) and put the zmatrix id numbers
c     of the two atoms in funcat1 and funcat2 Then simply include
c     this file where you want the energy to be calculated...
c
c-- If you don't understand this, then look at some examples

c***************************************************************
c       Work out the non-bonded parameters to use for this pair
c***************************************************************
c     interaction between constant clj parameters
      qiqj = ChgIJ(CLJPair(ipar1,ipar2))
      funcbij = BIJ(CLJPair(ipar1,ipar2))
      funcaij = AIJ(CLJPair(ipar1,ipar2))
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
      invfuncdist = ONE / dsqrt(funcdist2)
      invfuncdist6 = ONE / funcdist2**3
      invfuncdist12 = invfuncdist6**2


c***********************************************************
c       Now calculate the energies
c***********************************************************
c       **** reference state coulombic and LJ energies ****        
c        write(printstring,*) 'shift ',shift,' shiftf ',shiftf,' shiftb ',shiftb
c     .    ,' DeltaCore ',DeltaCore
c        call printLine(DEBUG,printstring)
      funccnrg = qiqj * invfuncdist
      funcljnrg = funcaij*invfuncdist12 - funcbij*invfuncdist6
c        write(6,*)"f(r)= ",ONE/invfuncdist6**(ONE/6.0d0), "g(r)= ",ONE/invfuncdist
c$$$      if (DifferentLambdas) then
c$$$c     *** Forwards/Backwards state coulombic and LJ energies *** 
c$$$         funccnrgf = qiqjf * invfuncdistf
c$$$         funcljnrgf = funcaijf*invfuncdist12f - funcbijf*invfuncdist6f
c$$$         funccnrgb = qiqjb * invfuncdistb
c$$$         funcljnrgb = funcaijb*invfuncdist12b - funcbijb*invfuncdist6b
c$$$c     write(6,*)"f(r)= ",ONE/invfuncdist6f**(ONE/6.0d0), "g(r)= ",ONE/invfuncdistf
c$$$c     write(6,*)"f(r)= ",ONE/invfuncdist6b**(ONE/6.0d0), "g(r)= ",ONE/invfuncdistb
c$$$      else
      funccnrgf = funccnrg
      funcljnrgf = funcljnrg
      funccnrgb = funccnrg
      funcljnrgb = funcljnrg 
c$$$      endif
      
c     funcmindist = funcdist2
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



    

