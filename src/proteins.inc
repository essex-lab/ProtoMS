
C  This file contains all of the variables that hold the protein structure
C  and parameters

c     The number of loaded proteins
      integer NProteins
c     The number of residues in each protein
      integer NResidues(MAXPROTEINS)
      
c     the name of the protein
      character*120 ProteinName(MAXPROTEINS)
c     the residue number of this residues
      integer ResNum(MAXPROTEINS,MAXRESIDUES)
c     the residue name of this residue
      character*4 ResNam(MAXPROTEINS,MAXRESIDUES)
c     the coordinates of the four backbone atoms
      double precision BBCoord(MAXPROTEINS,MAXRESIDUES,4,3)
c     the names of the four backbone atoms
      character*4 BBNam(MAXPROTEINS,MAXRESIDUES,4)
c     the atom number of the backbone atoms (the pdb atom number)
      integer BBNum(MAXPROTEINS,MAXRESIDUES,4)
c     the zmatrix atom number of the sidechain atoms...
      integer SCZid(MAXPROTEINS,MAXRESIDUES,MAXSCATOMS)
c     the same for the backbone atoms
      integer BBZid(MAXPROTEINS,MAXRESIDUES,4)
c     the coordinates of the sidechain atoms
      double precision SCCoord(MAXPROTEINS,MAXRESIDUES,MAXSCATOMS,3)
c     the atom name of one of these atoms
      character*4 SCNam(MAXPROTEINS,MAXRESIDUES,MAXSCATOMS)
c     the id number of the sidechain atoms (the atom number from the pdb file)
      integer SCNum(MAXPROTEINS,MAXRESIDUES,MAXSCATOMS)
c     the zmatrix atoms for the sidechain atoms
      integer ZMatSC(MAXPROTEINS,MAXRESIDUES,MAXSCATOMS,3)
c     the number of sidechain atoms for this protein and residue
      integer NSCAtoms(MAXPROTEINS,MAXRESIDUES)
c     the zmatrix values for the oxygen and sidechains
      double precision ZValO(MAXPROTEINS,MAXRESIDUES,3)
      double precision ZValSC(MAXPROTEINS,MAXRESIDUES,MAXSCATOMS,3)
c     whether or not these zmatrix parameters are variable
      logical ZMovO(MAXPROTEINS,MAXRESIDUES,3)
      logical ZMovSC(MAXPROTEINS,MAXRESIDUES,MAXSCATOMS,3)
c     whether or not to move a residue in the protein
      logical FixedRes(MAXPROTEINS,MAXRESIDUES)
c     whether or not to move the backbone of the residue in the protein
      logical FixedBac(MAXPROTEINS,MAXRESIDUES)
c     the acidic and basic residues in the protein
      integer AcidicRes(MAXPROTEINS,MAXACIDICRES)
      integer BasicRes(MAXPROTEINS,MAXBASICRES)
      integer NAcidicRes(MAXPROTEINS)
      integer NBasicRes(MAXPROTEINS)
c     the bonds between atoms (an index to the Bond array) - can include
c     inter-residue bonds, as use Bond array for whole molecule calculation,
c     so won't 'double count' them.
      integer ResidueBond(MAXPROTEINS,MAXRESIDUES,MAXBONDSPERRESIDUE)
c     the angle between atoms
      integer ResidueAngle(MAXPROTEINS,MAXRESIDUES,MAXANGLESPERRESIDUE)
c     the dihedrals between atoms
      integer ResidueDihedral(MAXPROTEINS,MAXRESIDUES,MAXDIHEDRALSPERRESIDUE)
c     the urey-bradley terms in each residue
      integer ResidueUreyBradley(MAXPROTEINS,MAXRESIDUES,MAXUREYBRADLEYSPERRESIDUE)
c     the number of bonds/angles/dihedrals in each residue...
      integer NResBonds(MAXPROTEINS,MAXRESIDUES)
      integer NResAngles(MAXPROTEINS,MAXRESIDUES)
      integer NResDihedrals(MAXPROTEINS,MAXRESIDUES)
      integer NResUreyBradleys(MAXPROTEINS,MAXRESIDUES)
c     the actual bonds between atoms (Zmat vals, + index to parameter
c      + residues involved)
      integer Bond(MAXPROTEINS,MAXBONDS,5)
      integer Angle(MAXPROTEINS,MAXANGLES,7)
      integer Dihedral(MAXPROTEINS,MAXDIHEDRALS,9)
      integer UreyBradley(MAXPROTEINS,MAXUREYBRADLEYS,7)
c     Whether or not these internals are dummy internals
      logical BondDummy(MAXPROTEINS,MAXBONDS)
      logical AngleDummy(MAXPROTEINS,MAXANGLES)
      logical DihedralDummy(MAXPROTEINS,MAXDIHEDRALS)
c     the list of residues that can move, and how many of them that are in
c     each protein
      integer ResMove(MAXPROTEINS,MAXRESIDUES)
      integer NResMove(MAXPROTEINS)
c     maximum amount to translate and rotate backbone of residue
      double precision ResADel(MAXPROTEINS,MAXRESIDUES)
      double precision ResRDel(MAXPROTEINS,MAXRESIDUES)
c     the number of bonded together residues in the protein
      integer NResResBonds(MAXPROTEINS)
c     the list of residues that are bonded together
      integer ResResBondIndex(MAXPROTEINS,MAXRESRESBONDS,2)
c     the inter and intra residue nb and 14 pair arrays
      integer NResNB(MAXPROTEINS,MAXRESIDUES)
      integer NResNB14(MAXPROTEINS,MAXRESIDUES)
      integer NResResNB(MAXPROTEINS,MAXRESRESBONDS)
      integer NResResNB14(MAXPROTEINS,MAXRESRESBONDS)
      integer ResNB(MAXPROTEINS,MAXRESIDUES,MAXRESNB,2)
      integer ResNB14(MAXPROTEINS,MAXRESIDUES,MAXRESNB14,2)
      integer ResResNB(MAXPROTEINS,MAXRESRESBONDS,MAXRESRESNB,2)
      integer ResResNB14(MAXPROTEINS,MAXRESRESBONDS,MAXRESRESNB14,2)
c     a working map of the residues...
      logical TempRes1Map(MAXSCATOMS+5,MAXSCATOMS+5)
      logical TempRes2Map(MAXSCATOMS+5,MAXSCATOMS+5)
      logical TempResResMap(MAXSCATOMS+5,MAXSCATOMS+5)
c     whether or not the bonds are flexible
      logical FlexBond(MAXPROTEINS,MAXBONDS)
      logical FlexAngle(MAXPROTEINS,MAXANGLES)
      logical FlexDihedral(MAXPROTEINS,MAXDIHEDRALS)
c     the delta for each bond/angle/dihedral
      double precision ZDeltaO(MAXPROTEINS,MAXRESIDUES,3)
      double precision ZDeltaSC(MAXPROTEINS,MAXRESIDUES,MAXSCATOMS,3)
c     the number of bonds/angles/dihedrals in the protein
      integer NBonds(MAXPROTEINS),NAngles(MAXPROTEINS),NDihedrals(MAXPROTEINS)
      integer NUreyBradleys(MAXPROTEINS)
c     the residue and atom number corresponding to the zmatrix atom number... (
c     sc from 1->MAXSCATOMS, BBAtoms from -1 to -4
      integer Zid2ResID(MAXPROTEINS,MAXATOMS,2)
c     a copy of all the coordinates indexed by protein and zmatrix id
      double precision ZCoords(MAXPROTEINS,MAXATOMS,3)
c     a copy of all the parameters indexde by protein and zmatrix id
      integer ZPar(MAXPROTEINS,MAXATOMS)
c     the total number of atoms in each protein
      integer NZAtoms(MAXPROTEINS)
c     the residues that are bonding together in the protein
      logical ResResBond(MAXPROTEINS,MAXRESIDUES,MAXRESIDUES)
c     the total charge on the protein
      double precision TotalCharge(MAXPROTEINS)
c     the geometric centre of each residue, and the radius**2
      double precision ResCentre(MAXPROTEINS,MAXRESIDUES,3)
      double precision ResRadius(MAXPROTEINS,MAXRESIDUES)
c     the geometric centre of each protein, and the old value for volume moves
      double precision ProCentre(MAXPROTEINS,3)
      double precision OldProCentre(MAXPROTEINS,3)
c     whether or not we have calculated the energy of each residue
      logical DoneResNrg(MAXPROTEINS,MAXRESIDUES)
c     the lengths/sizes of all of the bonds/angles/dihedrals...
      double precision BondLength(MAXPROTEINS,MAXBONDS)
      double precision AngleSize(MAXPROTEINS,MAXANGLES)
      double precision DihedralSize(MAXPROTEINS,MAXDIHEDRALS)
      double precision UreyBradleySize(MAXPROTEINS,MAXUREYBRADLEYS)
      logical ProteinScoop
c     whether or not to soften the residue-residue interactions
      logical ProteinSoftCore
      
      common /DPROTEINDATA/ BBCoord,SCCoord,ZValO,ZValSC,ZCoords,
     .                      ZDeltaO,ZDeltaSC,TotalCharge,
     .                      ResCentre,ResRadius,ProCentre,OldProCentre,
     .                      ResADel,ResRDel,BondLength,
     .                      AngleSize,DihedralSize,UreyBradleySize
      common /IPROTEINDATA/ NProteins,NResidues,ZMatSC,Zid2ResID,
     .                      NSCAtoms,ResNum,ResidueBond,BBNum,
     .                      ResidueAngle,ResidueDihedral,ResidueUreyBradley,
     .                      SCZid,BBZid,Bond,Angle,Dihedral,UreyBradley,
     .                      NResBonds,NResAngles,NResDihedrals,NResUreyBradleys,ZPar,
     .                      NBonds,NAngles,NDihedrals,NZAtoms,NUreyBradleys,
     .                      AcidicRes,BasicRes,NAcidicRes,NBasicRes,SCNum,
     .                      ResResBondIndex,NResResBonds,
     .                      NResNB,NResNB14,NResResNB,NResResNB14,
     .                      ResNB,ResNB14,ResResNB,ResResNB14,ResMove,NResMove
      common /LPROTEINDATA/ ZMovO,ZMovSC,FlexBond,FlexAngle,FlexDihedral,
     .                      BondDummy,AngleDummy,DihedralDummy,
     .                      ResResBond,TempRes1Map,TempRes2Map,TempResResMap,
     .                      FixedRes,FixedBac,DoneResNrg,ProteinScoop,
     .                      ProteinSoftCore
      common /CPROTEINDATA/ BBNam,SCNam,ResNam,ProteinName

      save /DPROTEINDATA/
      save /IPROTEINDATA/
      save /LPROTEINDATA/
      save /CPROTEINDATA/
      
