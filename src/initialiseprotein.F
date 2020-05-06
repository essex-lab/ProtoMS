      subroutine initialiseProtein(iprot)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'proteins.inc'
c####################################################
c
c     This code initialises the variables held 
c     for protein iprot - it effectly deletes 
c     protein iprot...
c
c     (C) Christopher Woods 24 January 2002
c
c####################################################

      integer iprot,i,j,k

      TotalCharge(iprot) = ZERO
      NAcidicRes(iprot) = 0
      NBasicRes(iprot) = 0
      NResResBonds(iprot) = 0
      ProteinName(iprot) = " "
      
      do i=1,MAXBASICRES
        BasicRes(iprot,i) = 0
      enddo
      
      do i=1,MAXACIDICRES
        AcidicRes(iprot,i) = 0
      enddo

      do i=1,MAXRESIDUES
        ResADel(iprot,i) = ZERO
        ResRDel(iprot,i) = ZERO
        NResBonds(iprot,i) = 0
        NResAngles(iprot,i) = 0
        NResDihedrals(iprot,i) = 0
        NResUreyBradleys(iprot,i) = 0
        ResNum(iprot,i) = i
        ResNam(iprot,i) = "   "
        NSCAtoms(iprot,i) = 0
        FixedRes(iprot,i) = .false.
        do j=1,MAXRESIDUES
          ResResBond(iprot,i,j) = .false.
        enddo
        do j=1,3
          do k=1,4
            BBCoord(iprot,i,k,j) = ZERO
          enddo
          do k=1,MAXSCATOMS
            SCCoord(iprot,i,k,j) = ZERO
            ZMatSC(iprot,i,k,j) = 0
            ZValSC(iprot,i,k,j) = ZERO
            ZMovSC(iprot,i,k,j) = .false.
            ZDeltaSC(iprot,i,k,j) = ZERO
          enddo
        enddo
        do j=1,MAXSCATOMS
          SCZid(iprot,i,j) = 0
          SCNam(iprot,i,j) = "   "
          SCNum(iprot,i,j) = 0
        enddo
        do j=1,4
          BBNam(iprot,i,j) = "   "
          BBNum(iprot,i,j) = 0
          BBZid(iprot,i,j) = 0
        enddo
        do j=1,MAXBONDSPERRESIDUE
          ResidueBond(iprot,i,j) = 0
        enddo
        do j=1,MAXANGLESPERRESIDUE
          ResidueAngle(iprot,i,j) = 0
        enddo
        do j=1,MAXDIHEDRALSPERRESIDUE
          ResidueDihedral(iprot,i,j) = 0
        enddo
        do j=1,MAXUREYBRADLEYSPERRESIDUE
          ResidueUreyBradley(iprot,i,j) = 0
        enddo
      enddo

      do i=1,MAXBONDS
        do j=1,5
          Bond(iprot,i,j) = 0
        enddo
        BondLength(iprot,i) = ZERO
        FlexBond(iprot,i) = .false.
      enddo
      do i=1,MAXANGLES
        do j=1,7
          Angle(iprot,i,j) = 0
        enddo
        AngleSize(iprot,i) = ZERO
        FlexAngle(iprot,i) = .false.
      enddo
      do i=1,MAXDIHEDRALS
        do j=1,9
          Dihedral(iprot,i,j) = 0
        enddo
        DihedralSize(iprot,i) = ZERO
        FlexDihedral(iprot,i) = .false.
      enddo
      do i=1,MAXUREYBRADLEYS
        do j=1,5
          UreyBradley(iprot,i,j) = 0
        enddo
        UreyBradleySize(iprot,i) = ZERO
      enddo

      do j=1,2
        do i=1,MAXATOMS
          Zid2ResID(iprot,i,j) = 0
        enddo
      enddo

      NBonds(iprot) = 0
      NAngles(iprot) = 0
      NDihedrals(iprot) = 0
      NUreyBradleys(iprot) = 0

      NZAtoms(iprot) = 0

      do i=1,MAXATOMS
        do j=1,3
          ZCoords(iprot,i,j) = ZERO
        enddo
        ZPar(iprot,i) = 0
      enddo

      end
