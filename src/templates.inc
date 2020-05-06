
c   This file holds the variables used for the templates

c     the names of atoms in the bonds/angles/dihedrals of the templates...
      character*4 TemplateBndName(MAXTEMPLATES,MAXTEMPLATEBONDS,2)
      character*4 TemplateAngName(MAXTEMPLATES,MAXTEMPLATEANGLES,3)
      character*4 TemplateDihName(MAXTEMPLATES,MAXTEMPLATEDIHEDRALS,4)
      character*4 TemplateUByName(MAXTEMPLATES,MAXTEMPLATEUREYBRADLEYS,3)
c     the actual residue name of the template
      character*4 TemplateName(MAXTEMPLATES)
      character*4 TemplateAtmName(MAXTEMPLATES,MAXTEMPLATEATMS,4)
c     the same for the molecule templates
      character*4 MolTemplateBndName(MAXMOLTEMPLATES,MAXMOLTEMPLATEBONDS,2)
      character*4 MolTemplateAngName(MAXMOLTEMPLATES,MAXMOLTEMPLATEANGLES,3)
      character*4 MolTemplateDihName(MAXMOLTEMPLATES,MAXMOLTEMPLATEDIHEDRALS,4)
      character*4 MolTemplateUByName(MAXMOLTEMPLATES,MAXMOLTEMPLATEUREYBRADLEYS,3)
      character*4 MolTemplateAtmRes(MAXMOLTEMPLATES,MAXMOLTEMPLATEATMS)
      character*4 MolTemplateBndRes(MAXMOLTEMPLATES,MAXMOLTEMPLATEBONDS,2)
      character*4 MolTemplateAngRes(MAXMOLTEMPLATES,MAXMOLTEMPLATEANGLES,3)
      character*4 MolTemplateDihRes(MAXMOLTEMPLATES,MAXMOLTEMPLATEDIHEDRALS,4)
      character*4 MolTemplateUByRes(MAXMOLTEMPLATES,MAXMOLTEMPLATEUREYBRADLEYS,3)
      character*4 MolTemplVarRes(MAXMOLTEMPLATES,MAXMOLTEMPLVAR)
      character*4 MolTemplateZRes(MAXMOLTEMPLATES,MAXMOLTEMPLATEATMS,3)
      character*20 MolTemplateName(MAXMOLTEMPLATES)
      character*4 MolTemplateAtmName(MAXMOLTEMPLATES,MAXMOLTEMPLATEATMS)
      character*4 MolTemplateZNam(MAXMOLTEMPLATES,MAXMOLTEMPLATEATMS,3)
      character*4 MolTemplVarNam(MAXMOLTEMPLATES,MAXMOLTEMPLVAR)
c     now the same for the solvent templates
      character*4 SvnTemplateName(MAXSVNTEMPLATES)
      character*4 GCSTemplateName(MAXSVNTEMPLATES)
      character*4 SvnTemplateAtmName(MAXSVNTEMPLATES,MAXSVNTEMPLATEATMS)
      character*4 GCSTemplateAtmName(MAXSVNTEMPLATES,MAXMOLTEMPLATEATMS)

      common /TEMPLATENAMES/ TemplateBndName,TemplateAngName,
     .     TemplateDihName,TemplateUByName,TemplateName,
     .     TemplateAtmName,MolTemplateBndName,
     .     MolTemplateAngName,MolTemplateDihName,MolTemplateUByName,MolTemplateName,
     .     MolTemplateAtmName,MolTemplateZNam,SvnTemplateName,GCSTemplateName,
     .     SvnTemplateAtmName,MolTemplVarNam,MolTemplateAtmRes,GCSTemplateAtmName,
     .     MolTemplateBndRes,MolTemplateAngRes,MolTemplateDihRes,MolTemplateUByRes,
     .     MolTemplateZRes,MolTemplVarRes

      save /TEMPLATENAMES/

c     template bond/angle/dihedral delta parameters
      double precision TemplateBondDelta(MAXTEMPLATES,MAXTEMPLATEBONDS)
      double precision TemplateAngleDelta(MAXTEMPLATES,MAXTEMPLATEANGLES)
      double precision TemplateDihedralDelta(MAXTEMPLATES,MAXTEMPLATEDIHEDRALS)
      double precision TemplateResADel(MAXTEMPLATES)
      double precision TemplateResRDel(MAXTEMPLATES)
c     the default zmatrix values for the residue and chain templates
      double precision TemplateZMatVal(MAXTEMPLATES,MAXTEMPLATEATMS,3)
      double precision ChainZMatVal(MAXCHAINTEMPLATES,MAXCHAINATOMS,3)
      double precision MolTemplateBondDelta(MAXMOLTEMPLATES,MAXMOLTEMPLATEBONDS)
      double precision MolTemplateAngleDelta(MAXMOLTEMPLATES,MAXMOLTEMPLATEANGLES)
      double precision MolTemplateDihedralDelta(MAXMOLTEMPLATES,MAXMOLTEMPLATEDIHEDRALS)
      double precision MolTemplVarVal(MAXMOLTEMPLATES,MAXMOLTEMPLVAR,2)
      double precision MolTemplateADel(MAXMOLTEMPLATES),MolTemplateRDel(MAXMOLTEMPLATES)
c     SG
      double precision MolTemplateTDel(MAXMOLTEMPLATES)
      double precision MolTemplateJCorr(MAXMOLTEMPLATES)
      integer MolTemplateJPmfOrd(MAXMOLTEMPLATES)
      integer MAXJPMFORDER
      parameter(MAXJPMFORDER=5)
      double precision MolTemplateJPmfCoeff(MAXMOLTEMPLATES,MAXJPMFORDER)
c     template solvent parameters
      double precision SvnTemplateADel(MAXSVNTEMPLATES)
      double precision SvnTemplateRDel(MAXSVNTEMPLATES)
      double precision GCSTemplateADel(MAXSVNTEMPLATES)
      double precision GCSTemplateRDel(MAXSVNTEMPLATES)
      integer          SvnTemplateOptimization(MAXSVNTEMPLATES)

      common /DTEMPLATEVALS/ TemplateAngleDelta,TemplateBondDelta,
     .           TemplateDihedralDelta,
     .           TemplateZmatVal,
     .           MolTemplateBondDelta,MolTemplateAngleDelta,
     .           MolTemplateDihedralDelta,SvnTemplateADel,GCSTemplateADel,GCSTemplateRDel,
     .           SvnTemplateRDel,MolTemplVarVal,MolTemplateAdel,MolTemplateRdel,
     .           MolTemplateTDel,MolTemplateJCorr,MolTemplateJPmfOrd,MolTemplateJPmfCoeff,
     .           TemplateResADel,TemplateResRDel,SvnTemplateOptimization

      save /DTEMPLATEVALS/

c     the number of bonds/angles/dihedrals in a particular template
      integer NTemplateBonds(MAXTEMPLATES)
      integer NTemplateAngles(MAXTEMPLATES)
      integer NTemplateDihedrals(MAXTEMPLATES)
      integer NTemplateUreyBradleys(MAXTEMPLATES)
      integer NTemplateAtms(MAXTEMPLATES)
      integer TemplateAtmCLJ(MAXTEMPLATES,MAXTEMPLATEATMS,6)
      integer TemplateBondPar(MAXTEMPLATES,MAXTEMPLATEBONDS,2)
      integer TemplateAnglePar(MAXTEMPLATES,MAXTEMPLATEANGLES,2)
      integer TemplateDihedralPar(MAXTEMPLATES,MAXTEMPLATEDIHEDRALS,2)
      integer TemplateUreyBradleyPar(MAXTEMPLATES,MAXTEMPLATEUREYBRADLEYS,2)
c     now the same for the molecule templates
      integer NMolTemplateBonds(MAXMOLTEMPLATES)
      integer NMolTemplateAngles(MAXMOLTEMPLATES)
      integer NMolTemplateDihedrals(MAXMOLTEMPLATES)
      integer NMolTemplateUreyBradleys(MAXMOLTEMPLATES)
      integer NMolTemplateAtms(MAXMOLTEMPLATES)
      integer MolTemplateBndPar(MAXMOLTEMPLATES,MAXMOLTEMPLATEBONDS,2)
      integer MolTemplateAngPar(MAXMOLTEMPLATES,MAXMOLTEMPLATEANGLES,2)
      integer MolTemplateDihPar(MAXMOLTEMPLATES,MAXMOLTEMPLATEDIHEDRALS,2)
      integer MolTemplateUByPar(MAXMOLTEMPlATES,MAXMOLTEMPLATEUREYBRADLEYS,2)
      integer MolTemplateAtmCLJ(MAXMOLTEMPLATES,MAXMOLTEMPLATEATMS,2)
c     JM Mod
      double precision MolTemplateDihWell(MAXMOLTEMPLATES,MAXMOLTEMPLATEDIHEDRALS,4)
      integer MolTemplateDihNWells(MAXMOLTEMPLATES,MAXMOLTEMPLATEDIHEDRALS)
      integer NMolTemplVar(MAXMOLTEMPLATES)
      integer MolTemplVarTyp(MAXMOLTEMPLATES,MAXMOLTEMPLVAR)
c     now the same for the solvent templates
      integer NSvnTemplateAtms(MAXSVNTEMPLATES)
      integer NGCSTemplateAtms(MAXSVNTEMPLATES)
      integer SvnTemplateAtmCLJ(MAXSVNTEMPLATES,MAXSVNTEMPLATEATMS,2)
      integer GCSTemplateAtmCLJ(MAXSVNTEMPLATES,MAXMOLTEMPLATEATMS,2)
c     the actual number of templates...
      integer NTemplates,NMolTemplates,NSvnTemplates,NGCSTemplates

      common /ITEMPLATEVALS/ NTemplateBonds,NTemplateAngles,
     .          NTemplateDihedrals,NTemplateUreyBradleys,
     .          NMolTemplateAtms,MolTemplateAtmCLJ,
     .          NMolTemplates,NMolTemplVar,MolTemplVarTyp,
     .          NTemplates,NTemplateAtms,TemplateAtmCLJ,
     .          NMolTemplateBonds,NMolTemplateAngles,NMolTemplateDihedrals,
     .          NMolTemplateUreyBradleys,MolTemplateUByPar,NGCSTemplateAtms,
     .          NSvnTemplateAtms,SvnTemplateAtmCLJ,NSvnTemplates,GCSTemplateAtmCLJ,
     .          MolTemplateBndPar,MolTemplateAngPar,MolTemplateDihPar,
     .          TemplateBondPar,TemplateAnglePar,TemplateDihedralPar,
     .          TemplateUreyBradleyPar,MolTemplateDihNWells,NGCSTemplates

      save /ITEMPLATEVALS/

      common /DMOLTEMPLATEVALS/ MolTemplateDihWell
      save /DMOLTEMPLATEVALS/
      
c     whether or not a particular bond/angle/dihedral is flexible and dummy
      logical TemplateBondType(MAXTEMPLATES,MAXTEMPLATEBONDS,2)
      logical TemplateAngleType(MAXTEMPLATES,MAXTEMPLATEANGLES,2)
      logical TemplateDihedralType(MAXTEMPLATES,MAXTEMPLATEDIHEDRALS,2)
      logical TemplateUreyBradleyType(MAXTEMPLATES,MAXTEMPLATEUREYBRADLEYS)
      logical MolTemplateBondType(MAXMOLTEMPLATES,MAXMOLTEMPLATEBONDS,2)
      logical MolTemplateAngleType(MAXMOLTEMPLATES,MAXMOLTEMPLATEANGLES,2)
      logical MolTemplateDihedralType(MAXMOLTEMPLATES,MAXMOLTEMPLATEDIHEDRALS,2)
      logical MolTemplateUreyBradleyType(MAXMOLTEMPLATES,MAXMOLTEMPLATEUREYBRADLEYS)
      
      common /LTEMPLATEVALS/ TemplateBondType,
     .          TemplateAngleType,TemplateDihedralType,TemplateUreyBradleyType,
     .          MolTemplateBondType,MolTemplateAngleType,MolTemplateDihedralType,
     .          MolTemplateUreyBradleyType

      save /LTEMPLATEVALS/

c     parameters for chain templates
      integer NChainTemplates
      integer NChainAtoms(MAXCHAINTEMPLATES)
      integer NChainBonds(MAXCHAINTEMPLATES)
      integer NChainAngles(MAXCHAINTEMPLATES)
      integer NChainDihedrals(MAXCHAINTEMPLATES)
      integer NChainUreyBradleys(MAXCHAINTEMPLATES)
      integer ChainBBPar(MAXCHAINTEMPLATES,4,2)
      character*4 ChainBBAtom(MAXCHAINTEMPLATES,4)
      character*20 ChainName(MAXCHAINTEMPLATES)
      character*20 TemplateResChain(MAXTEMPLATES,4)
      character*4 ChainAtm(MAXCHAINTEMPLATES,MAXCHAINATOMS,4)
      integer ChainAtmPar(MAXCHAINTEMPLATES,MAXCHAINATOMS,2)
      character*4 ChainBndName(MAXCHAINTEMPLATES,MAXCHAINBONDS,2)
      integer ChainBondPar(MAXCHAINTEMPLATES,MAXCHAINBONDS,2)
      logical ChainBondType(MAXCHAINTEMPLATES,MAXCHAINBONDS,2)
      double precision ChainBondDelta(MAXCHAINTEMPLATES,MAXCHAINBONDS)
      character*4 ChainAngName(MAXCHAINTEMPLATES,MAXCHAINANGLES,3)
      integer ChainAnglePar(MAXCHAINTEMPLATES,MAXCHAINANGLES,2)
      logical ChainAngleType(MAXCHAINTEMPLATES,MAXCHAINANGLES,2)
      double precision ChainAngleDelta(MAXCHAINTEMPLATES,MAXCHAINANGLES)
      character*4 ChainDihName(MAXCHAINTEMPLATES,MAXCHAINDIHEDRALS,4)
      integer ChainDihedralPar(MAXCHAINTEMPLATES,MAXCHAINDIHEDRALS,2)
      logical ChainDihedralType(MAXCHAINTEMPLATES,MAXCHAINDIHEDRALS,2)
      double precision ChainDihedralDelta(MAXCHAINTEMPLATES,MAXCHAINDIHEDRALS)
      character*4 ChainUByName(MAXCHAINTEMPLATES,MAXCHAINUREYBRADLEYS,3)
      integer ChainUreyBradleyPar(MAXCHAINTEMPLATES,MAXCHAINUREYBRADLEYS,2)
      logical ChainUreyBradleyType(MAXCHAINTEMPLATES,MAXCHAINUREYBRADLEYS)
      
      common /CHAINVARS/ ChainBondDelta,ChainAngleDelta,ChainDihedralDelta,ChainZmatVal,
     .                   NChainTemplates,NChainAtoms,NChainBonds,NChainAngles,
     .                   NChainDihedrals,NChainUreyBradleys,ChainBBPar,
     .                   ChainAtmPar,ChainBondPar,ChainAnglePar,ChainDihedralPar,
     .                   ChainUreyBradleyPar,ChainBBAtom,ChainName,ChainAtm,
     .                   ChainBndName,ChainAngName,ChainDihName,ChainUByName,
     .                   TemplateResChain,ChainBondType,ChainAngleType,ChainDihedralType,
     .                   ChainUreyBradleyType
     
      save /CHAINVARS/
      
c     The variables that are used for the overrides (to override the template parameters
c     thus allowing structure to be separated from forcefield)

c     the total number of loaded override templates
      integer NOverrides
c     the number of override parameters
      integer NOverrideParams(MAXOVERRIDES)
c     the atom name of each override parameter
      character*4 OverrideAtm(MAXOVERRIDES,MAXOVERRIDEPARAMS)
c     the override parameters
      integer OverrideParam(MAXOVERRIDES,MAXOVERRIDEPARAMS,2)

c     the index into the override template array for the overrides for each residue
c     in a particular chain position (5 possible positions)
      integer ResOverrideIndex(5,MAXTEMPLATES)
c     the total number of override templates for each chain position (5 available positions)
      integer NResOverrides(5)
c     the name of the residue in each chain position
      character*4 ResOverrideNam(5,MAXTEMPLATES)
                  
      common /OVERRIDEVARS/ NOverrides,NOverrideParams,OverrideParam,ResOverrideIndex,
     .                      NResOverrides,ResOverrideNam,OverrideAtm
      
      save /OVERRIDEVARS/
      
