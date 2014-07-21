      subroutine initialiseSolute(isol)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'solutes.inc'
      include 'moves.inc'
c####################################################
c
c     This code initialises the variables held 
c     for solute isol - it effectly deletes 
c     solute isol...
c
c     (C) Christopher Woods 28 January 2002
c
c####################################################

      integer isol,i,j

      SoluteName(isol) = " "
      NSolAtoms(isol) = 0
      do i=1,MAXSOLUTEATOMS
        do j=1,3
          ZSolCoords(isol,i,j) = ZERO
          ZSolMat(isol,i,j) = 0
          ZSolVal(isol,i,j) = ZERO
          ZSolDelta(isol,i,j) = ZERO
          ZSolMov(isol,i,j) = .false.
        enddo
        ZSolPar(isol,i) = 0
        SolAtomNam(isol,i) = NULLNAME
        SolAtomNum(isol,i) = 0
      enddo

      DualTopologySolute(isol) = 0
      
      do i=1,MAXSOLUTERESIDUES
        SolResidues(isol,i,1) = 0
        SolResidues(isol,i,2) = 0
        SolResNam(isol,i) = NULLNAME
        SolResNum(isol,i) = 0
        NSolTry(isol,i) = 0
        NSolAcp(isol,i) = 0
      enddo
      NSolResidues(isol) = 0
      
      do i=1,MAXSOLUTEBONDS
        do j=1,3
          SolBond(isol,i,j) = 0 
        enddo
        SolBondLength(isol,i) = ZERO
        SolBondDelta(isol,i) = ZERO
        FlexSolBond(isol,i) = .false.
        SolBondDummy(isol,i) = .false.
      enddo
      NSolBonds(isol) = 0

      do i=1,MAXSOLUTEANGLES
        do j=1,4
          SolAngle(isol,i,j) = 0 
        enddo
        SolAngleSize(isol,i) = ZERO
        SolAngleDelta(isol,i) = ZERO
        FlexSolAngle(isol,i) = .false.
        SolAngleDummy(isol,i) = .false.
      enddo
      NSolAngles(isol) = 0

      do i=1,MAXSOLUTEDIHEDRALS
        do j=1,5
          SolDihedral(isol,i,j) = 0 
        enddo
        SolDihedralSize(isol,i) = ZERO
        SolDihedralDelta(isol,i) = ZERO
        FlexSolDihedral(isol,i) = .false.
        SolDihedralDummy(isol,i) = .false.
      enddo
      NSolDihedrals(isol) = 0

      do i=1,MAXSOLUTEUREYBRADLEYS
        do j=1,3
          SolUreyBradley(isol,i,j) = 0
        enddo
        SolUreyBradleySize(isol,i) = ZERO
      enddo
      NSolUreyBradleys(isol) = 0

      NSolNB(isol) = 0
      do i=1,MAXSOLUTENBPAIRS
        SolNB(isol,i,1) = 0
        SolNB(isol,i,2) = 0
        SolNB14(isol,i) = .false.
      enddo
      
      end
