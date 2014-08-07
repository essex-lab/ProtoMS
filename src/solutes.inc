
c   This file has all of the variables that hold the solute coordinates and parameters

c     the number of loaded solutes
      integer NSolutes

c     the names of the solutes...
      character*120 SoluteName(MAXSOLUTES)

c     the number of solute groups
      integer NSolGroups
c     the number of solutes in each group
      integer NGroupSols(MAXSOLUTES)
c     the list of solutes in each group
      integer SolGroup(MAXSOLUTES,MAXSOLUTES)
c     the number of atoms in the solute
      integer NSolAtoms(MAXSOLUTES)
c     cartesian coordinates
      double precision ZSolCoords(MAXSOLUTES,MAXSOLUTEATOMS,3)
c     the coordinates of the three dummy atoms
      double precision SolDummy(MAXSOLUTES,3,3)
c     the clj parameter types
      integer ZSolPar(MAXSOLUTES,MAXSOLUTEATOMS)
c     the zmatrix atoms
      integer ZSolMat(MAXSOLUTES,MAXSOLUTEATOMS,3)
c     whether or not to move this zmatrix coordinate
      logical ZSolMov(MAXSOLUTES,MAXSOLUTEATOMS,3)
c     the delta value with which to move the coordinates
      double precision ZSolDelta(MAXSOLUTES,MAXSOLUTEATOMS,3)
c     the zmatrix values
      double precision ZSolVal(MAXSOLUTES,MAXSOLUTEATOMS,3)
c     JM mod, wells for some dihedrals..
      integer ZSolDihNWells(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision ZSolDihWell(MAXSOLUTES,MAXSOLUTEATOMS,4)
c     the start and end indicies of the residues...
      integer NSolResidues(MAXSOLUTES)
      integer SolResidues(MAXSOLUTES,MAXSOLUTERESIDUES,2)
      integer SolResidue(MAXSOLUTES,MAXSOLUTEATOMS)
c     the bonds/angles/dihedrals in the solutes
      integer SolBond(MAXSOLUTES,MAXSOLUTEBONDS,3)
      integer GCSBond(125,MAXSOLUTEBONDS,3)
      integer SolAngle(MAXSOLUTES,MAXSOLUTEANGLES,4)
      integer GCSAngle(125,MAXSOLUTEANGLES,4)
      integer SolDihedral(MAXSOLUTES,MAXSOLUTEDIHEDRALS,5)
      integer GCSDihedral(125,MAXSOLUTEDIHEDRALS,5)
      integer SolUreyBradley(MAXSOLUTES,MAXSOLUTEUREYBRADLEYS,4)
      integer GCSUreyBradley(125,MAXSOLUTEUREYBRADLEYS,4)
      integer SolNB(MAXSOLUTES,MAXSOLUTENBPAIRS,2)
      integer GCSNB(125,MAXSOLUTENBPAIRS,2)
      logical SolNB14(MAXSOLUTES,MAXSOLUTENBPAIRS)
      double precision SolBondDelta(MAXSOLUTES,MAXSOLUTEBONDS)
      double precision SolAngleDelta(MAXSOLUTES,MAXSOLUTEANGLES)
      double precision SolDihedralDelta(MAXSOLUTES,MAXSOLUTEDIHEDRALS)
      double precision SolBondLength(MAXSOLUTES,MAXSOLUTEBONDS)
      double precision GCSBondLength(125,MAXSOLUTEBONDS)
      double precision SolAngleSize(MAXSOLUTES,MAXSOLUTEANGLES)
      double precision GCSAngleSize(125,MAXSOLUTEANGLES)
      double precision SolDihedralSize(MAXSOLUTES,MAXSOLUTEDIHEDRALS)
      double precision GCSDihedralSize(125,MAXSOLUTEDIHEDRALS)
      double precision SolUreyBradleySize(MAXSOLUTES,MAXSOLUTEUREYBRADLEYS)
      double precision GCSUreyBradleySize(125,MAXSOLUTEUREYBRADLEYS)
      double precision SolADel(MAXSOLUTES),SolRDel(MAXSOLUTES)
c     SG
      double precision SolTDel(MAXSOLUTES)
      double precision SolJCorr(MAXMOLTEMPLATES)
      integer SolJPmfOrd(MAXMOLTEMPLATES)
      integer MAXORDER
      parameter(MAXORDER=5)
      double precision SolJPmfCoeff(MAXMOLTEMPLATES,MAXORDER)
      
      logical SolBondDummy(MAXSOLUTES,MAXSOLUTEBONDS)
      logical SolAngleDummy(MAXSOLUTES,MAXSOLUTEANGLES)
      logical SolDihedralDummy(MAXSOLUTES,MAXSOLUTEDIHEDRALS)
      integer NSolNB(MAXSOLUTES)
      integer NGCSNB(125)
      integer NSolBonds(MAXSOLUTES)
      integer NGCSBonds(125)
      integer NSolAngles(MAXSOLUTES)
      integer NGCSAngles(125)
      integer NSolDihedrals(MAXSOLUTES)
      integer NGCSDihedrals(125)
      integer NSolUreyBradleys(MAXSOLUTES)
      integer NGCSUreyBradleys(125)
      integer NSolVariable(MAXSOLUTES)
      integer NSolForwards(MAXSOLUTES)
      integer NSolBackwards(MAXSOLUTES)
      integer DualTopologySolute(MAXSOLUTES)
      logical DualTopologySyncTrans(MAXSOLUTES)
      logical DualTopologySyncRot(MAXSOLUTES)
      logical DualTopologySyncInt(MAXSOLUTES)
      double precision SolVariableVal(MAXSOLUTES,MAXSOLVARIABLES,2)
      integer SolVariableTyp(MAXSOLUTES,MAXSOLVARIABLES)
      integer SolVariableAtm(MAXSOLUTES,MAXSOLVARIABLES)
      logical IgnoreSolute(MAXSOLUTES)
c     whether or not the bonds/angles/dihedrals are flexible
      logical FlexSolBond(MAXSOLUTES,MAXSOLUTEBONDS)
      logical FlexSolAngle(MAXSOLUTES,MAXSOLUTEANGLES)
      logical FlexSolDihedral(MAXSOLUTES,MAXSOLUTEDIHEDRALS)
c     the names of the atoms...
      character*4 SolResNam(MAXSOLUTES,MAXSOLUTERESIDUES)
      character*4 SolAtomNam(MAXSOLUTES,MAXSOLUTEATOMS)
c     the id numbers of the solutes atoms and residues in the pdb file
      integer SolResNum(MAXSOLUTES,MAXSOLUTERESIDUES)
      integer SolAtomNum(MAXSOLUTES,MAXSOLUTEATOMS)
c     the solute centre of the solute, and the radius**2
      double precision SolCentre(MAXSOLUTES,3)
c     the old solute centre - used for cutoff
      double precision OldSolCentre(MAXSOLUTES,3)
      double precision SolRadius(MAXSOLUTES)
c     the solute geometric centre -- for rotation/translations
      double precision SolGeomCentre(MAXSOLUTES,3)
      double precision OldSolGeomCentre(MAXSOLUTES,3)
c     whether or not we have calculated the energy of this solute
      logical DoneSolNrg(MAXSOLUTES)
c     whether or not a softcore is used for interactions involving that solute
      logical SoluteSoftCore(MAXSOLUTES)
c     wheter to use library charges for intra molecular energies or the perturbed ones
      logical perturbeIntra

c     SG
c     scaling parameter for solute
      double precision SoluteTheta(MAXSOLUTES)  
      double precision SoluteHydScale(MAXSOLUTES)   
      logical UseConcCorr

c     Number of theta solute kinds
      integer NThetaSolKind
c     Number of theta solute molecules
      integer NThetaSolutes
c     Which kind of theta solute a solute is
      integer ThetaKind(MAXSOLUTES)

      common /DSOLDATA/ ZSolCoords,ZSolVal,SolBondLength,SolAngleSize,
     .                  SolDihedralSize,SolUreyBradleySize,SolBondDelta,SolAngleDelta,
     .                  SolDihedralDelta,ZSolDelta,SolDummy,
     .                  SolCentre,OldSolCentre,SolRadius,SolGeomCentre,
     .                  OldSolGeomCentre,SolADel,SolRDel,SolTDel,SolJCorr,SolJPmfOrd,SolJPmfCoeff,
     .                  SolVariableVal,ZSolDihWell
      common /ISOLDATA/ NSolutes,ZSolPar,ZSolMat,SolResidues,SolBond,SolAngle,
     .                  SolDihedral,SolUreyBradley,NSolBonds,NSolAngles,NSolDihedrals,
     .                  NSolUreyBradleys,NSolAtoms,NSolResidues,SolResidue,
     .                  SolNB,NSolNB,NSolGroups,NGroupSols,SolGroup,
     .                  SolVariableTyp,SolVariableAtm,NSolForwards,
     .                  NSolBackwards,NSolVariable,DualTopologySolute,
     .                  SolResNum,SolAtomNum,ZSolDihNWells
      common /GCSSDATA/ GCSBond,GCSAngle,GCSDihedral,GCSUreyBradley,GCSBondLength,NGCSBonds,
     .                  NGCSAngles,GCSNB,NGCSDihedrals,NGCSUreyBradleys,GCSAngleSize,GCSDihedralSize,
     .                  GCSUreyBradleySize,NGCSNB
      common /CSOLDATA/ SolResNam,SolAtomNam,SoluteName
      common /LSOLDATA/ FlexSolBond,FlexSolAngle,FlexSolDihedral,ZSolMov,DoneSolNrg,
     .                  SolNB14,IgnoreSolute,DualTopologySyncTrans,SoluteTheta,SoluteHydScale, 
     .                  DualTopologySyncRot,DualTopologySyncInt,
     .                  SolBondDummy,SolAngleDummy,SolDihedralDummy,SoluteSoftCore,perturbeIntra                  

      common /FRAGDATA/ NThetaSolKind,NThetaSolutes,ThetaKind,UseConcCorr

      save /DSOLDATA/
      save /ISOLDATA/
      save /CSOLDATA/
      save /LSOLDATA/
      save /GCSSDATA/
      save /FRAGDATA/
