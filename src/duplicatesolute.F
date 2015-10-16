      integer function duplicateSolute(isol)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'solutes.inc'
c#########################################################
c
c     This creates a copy of solute isol, and returns
c     the index of the new solute
c
c     (C) Christopher Woods   14 August 2002
c
c#########################################################
 
      integer isol,i,j
                
      duplicateSolute = 0
      
      if (isol.lt.0 .or. isol.gt.NSolutes) return

      NSolutes = NSolutes + 1
      if (NSolutes.gt.MAXSOLUTES) then
        call printLine(FATAL,"Exceeded maximum number of solutes - cannot duplicate!")
        write(printstring,*)"Increase MAXSOLUTES and recompile"
        call closeProgram(-1,printstring)
        return
      endif
      
c     copy the solute's name
      SoluteName(NSolutes) = SoluteName(isol)      
      
c     copy the bonds/angles/dihedrals
      NSolAtoms(NSolutes) = NSolAtoms(isol)
      NSolResidues(NSolutes) = NSolResidues(isol)
      do j=1,3
        do i=1,NSolAtoms(isol)           
          ZSolCoords(NSolutes,i,j) = ZSolCoords(isol,i,j)
          ZSolMat(NSolutes,i,j) = ZSolMat(isol,i,j)
          ZSolMov(NSolutes,i,j) = ZSolMov(isol,i,j)
          ZSolDelta(NSolutes,i,j) = ZSolDelta(isol,i,j)
          ZSolVal(NSolutes,i,j) = ZSolVal(isol,i,j)
        enddo
        do i=1,3
          SolDummy(NSolutes,i,j) = SolDummy(isol,i,j)
        enddo
      enddo

      do i=1,NSolAtoms(isol)
        ZSolPar(NSolutes,i) = ZSolPar(isol,i)
        SolResidue(NSolutes,i) = SolResidue(isol,i)
        SolAtomNam(NSolutes,i) = SolAtomNam(isol,i)
      enddo

      do i=1,NSolResidues(isol)
        SolResidues(NSolutes,i,1) = SolResidues(isol,i,1)
        SolResidues(NSolutes,i,2) = SolResidues(isol,i,2)
        SolResNam(NSolutes,i) = SolResNam(isol,i)
      enddo

      do i=1,3
        SolCentre(NSolutes,i) = SolCentre(isol,i)
      enddo
      SolRadius(NSolutes) = SolRadius(isol)

      NSolBonds(NSolutes) = NSolBonds(isol)
      NSolAngles(NSolutes) = NSolAngles(isol)
      NSolDihedrals(NSolutes) = NSolDihedrals(isol)
      NSolNB(NSolutes) = NSolNB(isol)
                      
      do i=1,NSolBonds(isol)
        do j=1,3
          SolBond(NSolutes,i,j) = SolBond(isol,i,j)
        enddo
        SolBondDelta(NSolutes,i) = SolBondDelta(isol,i)
        SolBondLength(NSolutes,i) = SolBondLength(isol,i)
        FlexSolBond(NSolutes,i) = FlexSolBond(isol,i)
      enddo

      do i=1,NSolAngles(isol)
        do j=1,4
          SolAngle(NSolutes,i,j) = SolAngle(isol,i,j)
        enddo
        SolAngleDelta(NSolutes,i) = SolAngleDelta(isol,i)
        SolAngleSize(NSolutes,i) = SolAngleSize(isol,i)
        FlexSolAngle(NSolutes,i) = FlexSolAngle(isol,i)
      enddo

      do i=1,NSolDihedrals(isol)
        do j=1,5
          SolDihedral(NSolutes,i,j) = SolDihedral(isol,i,j)
        enddo
        SolDihedralDelta(NSolutes,i) = SolDihedralDelta(isol,i)
        SolDihedralSize(NSolutes,i) = SolDihedralDelta(isol,i)
        FlexSolDihedral(NSolutes,i) = FlexSolDihedral(isol,i)
      enddo
      
      do i=1,NSolUreyBradleys(isol)
        do j=1,3
          SolUreyBradley(NSolutes,i,j) = SolUreyBradley(isol,i,j)
        enddo
        SolUreyBradleySize(NSolutes,i) = SolUreyBradleySize(isol,i)
      enddo
      
      do i=1,NSolNB(isol)
        do j=1,2
          SolNB(NSolutes,i,j) = SolNB(isol,i,j)
        enddo
        SolNB14(NSolutes,i) = SolNB14(isol,i)
      enddo

      duplicateSolute = NSolutes
              
      return
      end
