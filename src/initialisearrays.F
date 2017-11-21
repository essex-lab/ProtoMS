
      subroutine initialiseArrays
      include 'dimensions.inc'
      include 'metropolis.inc'
      include 'averages.inc'
      include 'flags.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'templates.inc'
      include 'constants.inc'
      include 'solvents.inc'
      include 'simparams.inc'
      include 'enums.inc'
      include 'moves.inc'
      include 'parameters.inc'
      include 'pressure.inc'
      include 'energies.inc'
      include 'stacks.inc'
      include 'printer.inc'
      include 'extraenergies.inc'
      include 'labels.inc'
      include 'rotamer.inc'
      include 'shrakeinit.inc'
      include 'quickpot.inc'
c      include 'cavity.inc'
      
      integer i,j
      double precision deg2Rad

c      integer k

c     set the default 1.0 / 1-4 scale factor
      OneOverScl14C = 0.5
      OneOverScl14LJ = 0.5

c     default is to use geometric combining rules
      lorentzberthelot = .FALSE.

c     set the MC bias to ONE (no bias)
      mt_bias = ONE

      NProteins = 0
      NSolutes = 0
      NSolvents = 0
      NGCSolutes = 0
      
c     set the big move probability for rotamer moves to 10%
      BigMoveProb = 0.1D+00

c     set the delta for big moves to +/- 180 degrees
      BigMoveDelta = deg2Rad(ONEEIGHTY)

c     set the delta for rotamer moves to +/- 10 degrees
      RotamerDelta = deg2Rad(10.0D+00)

      startstep = 0
      completed = .true.
      
      do i=1,3
        BoxDimension(i) = ZERO
        BoxDimension2(i) = ZERO
      enddo

      do i=1,MAXPROTEINS
        NResidues(i) = 0
        do j=1,3
          ProCentre(i,j) = ZERO
          OldProCentre(i,j) = ZERO
        enddo
        do j=1,MAXRESIDUES
          DoneResNrg(i,j) = .false. 
          NResTry(i,j) = 0
          NResAcp(i,j) = 0
        enddo
      enddo

      do i=1,MAXCLJTYPES
        do j=1,MAXCLJTYPES
          CLJPair(i,j) = 0
        enddo
        CLJType(i) = 0
      enddo
                
      do i=1,MAXCLJ
        GotCLJ(i,1) = .false.
        GotCLJ(i,2) = .false.
        ChgDB(i) = ZERO
        SigDB(i) = ZERO
        EpsDB(i) = ZERO
        AtmDB(i) = NULLCLJ
      enddo

      do i=1,MAXCLJPAIRS
        ChgIJ(i) = ZERO
        AIJ(i) = ZERO
        BIJ(i) = ZERO
      enddo

      NeutraliseProtein = .false.

      NCLJTypes = 0

      do i=1,MAXBNDPARAM
        do j=1,2
          BndParam(i,j) = ZERO
          BndName(i,j) = NULLNAME
        enddo
        GotBndParam(i) = .false.
        BndParamIndex(i) = 0
      enddo

      do i=1,MAXANGPARAM
        do j=1,2
          AngParam(i,j) = ZERO
          AngName(i,j) = NULLNAME
        enddo
        GotAngParam(i) = .false.
        AngName(i,3) = NULLNAME
        AngParamIndex(i) = 0
      enddo

c     the first bond/angle/dihedral parameter is the null parameter...
      NAngParam = 0
      NBndParam = 0
      NDihParam = 0
      NUByParam = 0
    
      do i=1,MAXDIHPARAM
        do j=1,4
          DihName(i,j) = NULLNAME
        enddo
        GotDihParam(i) = .false.
        NDihTerms(i) = 0
        do j=1,MAXDIHTERMSPERDIHEDRAL
          DihTermList(i,j) = 0
        enddo
        DihParamIndex(i) = 0
      enddo
      
      do i=1,MAXDIHTERMS
        do j=1,4
          DihTerm(i,j) = ZERO
        enddo
        GotDihTerm(i) = .false.
      enddo
      
      do i=1,MAXUBYPARAM
        do j=1,2
          UByParam(i,j) = ZERO
          UByName(i,j) = NULLNAME
        enddo
        GotUByParam(i) = .false.
        UByParamIndex(i) = 0
      enddo
        
c     initialising the templates..
      do i=1,MAXTEMPLATES
        do j=1,MAXTEMPLATEBONDS
          TemplateBndName(i,j,1) = NULLNAME
          TemplateBndName(i,j,2) = NULLNAME
          TemplateBondType(i,j,1) = .false.
          TemplateBondType(i,j,2) = .false.
          TemplateBondDelta(i,j) = ZERO
          TemplateBondPar(i,j,1) = -1
          TemplateBondPar(i,j,2) = -1
        enddo
        NTemplateBonds(i) = 0
        do j=1,MAXTEMPLATEANGLES
          TemplateAngName(i,j,1) = NULLNAME
          TemplateAngName(i,j,2) = NULLNAME
          TemplateAngName(i,j,3) = NULLNAME
          TemplateAngleType(i,j,1) = .false.
          TemplateAngleType(i,j,2) = .false.
          TemplateAngleDelta(i,j) = ZERO
          TemplateAnglePar(i,j,1) = -1
          TemplateAnglePar(i,j,2) = -1
        enddo
        TemplateResADel(i) = ZERO
        TemplateResRDel(i) = ZERO
        NTemplateAngles(i) = 0
        do j=1,MAXTEMPLATEDIHEDRALS
          TemplateDihName(i,j,1) = NULLNAME          
          TemplateDihName(i,j,2) = NULLNAME          
          TemplateDihName(i,j,3) = NULLNAME          
          TemplateDihName(i,j,4) = NULLNAME          
          TemplateDihedralType(i,j,1) = .false.
          TemplateDihedralType(i,j,2) = .false.
          TemplateDihedralDelta(i,j) = ZERO
          TemplateDihedralPar(i,j,1) = -1
          TemplateDihedralPar(i,j,2) = -1
        enddo
        NTemplateDihedrals(i) = 0
        do j=1,MAXTEMPLATEUREYBRADLEYS
          TemplateUByName(i,j,1) = NULLNAME
          TemplateUByName(i,j,2) = NULLNAME
          TemplateUreyBradleyType(i,j) = .false.
          TemplateUreyBradleyPar(i,j,1) = -1
          TemplateUreyBradleyPar(i,j,2) = -1
        enddo
        NTemplateUreyBradleys(i) = 0
        do j=1,MAXTEMPLATEATMS
          TemplateAtmName(i,j,1) = NULLNAME
          TemplateAtmName(i,j,2) = NULLNAME
          TemplateAtmName(i,j,3) = NULLNAME
          TemplateAtmName(i,j,4) = NULLNAME
          TemplateAtmCLJ(i,j,1) = 0
          TemplateAtmCLJ(i,j,2) = 0
        enddo
        NTemplateAtms(i) = 0
        NTemplates = 0
      enddo

      NChainTemplates = 0
      do i=1,MAXCHAINTEMPLATES
        NChainAtoms(i) = 0
        NChainBonds(i) = 0
        NChainAngles(i) = 0
        NChainDihedrals(i) = 0
        NChainUreyBradleys(i) = 0
        do j=1,4
          ChainBBPar(i,j,1) = 0
          ChainBBPar(i,j,2) = 0
          ChainBBAtom(i,j) = '   '
        enddo
      enddo                    
c
c     This is a temporary hack to get gcc to compile 
c      
c$$$      do i=1,MAXMOLTEMPLATES
c$$$        do j=1,MAXMOLTEMPLATEBONDS
c$$$          MolTemplateBndName(i,j,1) = NULLNAME
c$$$          MolTemplateBndName(i,j,2) = NULLNAME
c$$$          MolTemplateBndRes(i,j,1) = NULLNAME
c$$$          MolTemplateBndRes(i,j,2) = NULLNAME
c$$$          MolTemplateBondType(i,j,1) = .false.
c$$$          MolTemplateBondType(i,j,2) = .false.
c$$$          MolTemplateBondDelta(i,j) = ZERO
c$$$          MolTemplateBndPar(i,j,1) = -1
c$$$          MolTemplateBndPar(i,j,2) = -1
c$$$        enddo
c$$$        NMolTemplateBonds(i) = 0
c$$$        do j=1,MAXMOLTEMPLATEANGLES
c$$$          MolTemplateAngName(i,j,1) = NULLNAME
c$$$          MolTemplateAngName(i,j,2) = NULLNAME
c$$$          MolTemplateAngName(i,j,3) = NULLNAME
c$$$          MolTemplateAngRes(i,j,1) = NULLNAME
c$$$          MolTemplateAngRes(i,j,2) = NULLNAME
c$$$          MolTemplateAngRes(i,j,3) = NULLNAME
c$$$          MolTemplateAngleType(i,j,1) = .false.
c$$$          MolTemplateAngleType(i,j,2) = .false.
c$$$          MolTemplateAngleDelta(i,j) = ZERO
c$$$          MolTemplateAngPar(i,j,1) = -1
c$$$          MolTemplateAngPar(i,j,2) = -1
c$$$        enddo
c$$$        NMolTemplateAngles(i) = 0
c$$$        do j=1,MAXMOLTEMPLATEDIHEDRALS
c$$$          MolTemplateDihName(i,j,1) = NULLNAME
c$$$          MolTemplateDihName(i,j,2) = NULLNAME
c$$$          MolTemplateDihName(i,j,3) = NULLNAME
c$$$          MolTemplateDihName(i,j,4) = NULLNAME
c$$$          MolTemplateDihRes(i,j,1) = NULLNAME
c$$$          MolTemplateDihRes(i,j,2) = NULLNAME
c$$$          MolTemplateDihRes(i,j,3) = NULLNAME
c$$$          MolTemplateDihRes(i,j,4) = NULLNAME
c$$$          MolTemplateDihedralType(i,j,1) = .false.
c$$$          MolTemplateDihedralType(i,j,2) = .false.
c$$$          MolTemplateDihedralDelta(i,j) = ZERO
c$$$          MolTemplateDihPar(i,j,1) = -1
c$$$          MolTemplateDihPar(i,j,2) = -1
c$$$        enddo
c$$$        NMolTemplateDihedrals(i) = 0
c$$$        do j=1,MAXMOLTEMPLATEUREYBRADLEYS
c$$$          MolTemplateUByName(i,j,1) = NULLNAME
c$$$          MolTemplateUByName(i,j,2) = NULLNAME
c$$$          MolTemplateUBYName(i,j,3) = NULLNAME
c$$$          MolTemplateUByRes(i,j,1) = NULLNAME
c$$$          MolTemplateUByRes(i,j,2) = NULLNAME
c$$$          MolTemplateUByRes(i,j,3) = NULLNAME
c$$$          MolTemplateUreyBradleyType(i,j) = .false.
c$$$          MolTemplateUByPar(i,j,1) = -1
c$$$          MolTemplateUByPar(i,j,2) = -1
c$$$        enddo
c$$$        NMolTemplateUreyBradleys(i) = 0                    
c$$$        do j=1,MAXMOLTEMPLATEATMS
c$$$          MolTemplateAtmName(i,j) = NULLNAME
c$$$          MolTemplateAtmRes(i,j) = NULLNAME
c$$$          MolTemplateAtmCLJ(i,j,1) = 0
c$$$          MolTemplateAtmCLJ(i,j,2) = 0
c$$$        enddo
c$$$        NMolTemplateAtms(i) = 0   
c$$$        NMolTemplVar(i) = 0
c$$$        do j=1,MAXMOLTEMPLVAR
c$$$          MolTemplVarNam(i,j) = NULLNAME
c$$$          MolTemplVarRes(i,j) = NULLNAME
c$$$          MolTemplVarTyp(i,j) = 0
c$$$          MolTemplVarVal(i,j,1) = ZERO
c$$$          MolTemplVarVal(i,j,2) = ZERO
c$$$        enddo
c$$$        MolTemplateADel(i) = ZERO
c$$$        MolTemplateRDel(i) = ZERO
c$$$      enddo
      NMolTemplates = 0
      do i=1,MAXSVNTEMPLATES
        SvnTemplateADel(i) = 0.0
        SvnTemplateRDel(i) = 0.0
        SvnTemplateName(i) = NULLNAME
        NSvnTemplateAtms(i) = 0
        NGCSTemplateAtms(i) = 0
        do j=1,MAXSVNTEMPLATEATMS
          SvnTemplateAtmName(i,j) = NULLNAME
          SvnTemplateAtmCLJ(i,j,1) = 0
          SvnTemplateAtmCLJ(i,j,2) = 0
        enddo
      enddo
      NSvnTemplates = 0
      NGCSTemplates = 0
c      NGCSTemplateAtmName(MAXSVNTEMPLATES,MAXMOLTEMPLATES) = 0
c     initialise the energies...    
      NewTotalEnergy = ZERO
      OldTotalEnergy = ZERO
      NewSvnSvnC = ZERO
      NewSvnSvnLJ = ZERO
      NewSvnSvnDist = ZERO
      OldSvnSvnC = ZERO
      OldSvnSvnLJ = ZERO
      OldSvnSvnDist = ZERO
      AvgSvnSvnC = ZERO
      AvgSvnSvnC2 = ZERO
      AvgSvnSvnLJ = ZERO
      AvgSvnSvnLJ2 = ZERO
      
      NewGCSGCSC = ZERO
      NewGCSGCSLJ = ZERO
      OldGCSGCSC = ZERO
      OldGCSGCSLJ = ZERO
      AvgGCSGCSC = ZERO
      AvgGCSGCSC2 = ZERO
      AvgGCSGCSLJ = ZERO
      AvgGCSGCSLJ2 = ZERO

      do i=1,MAXPROTEINS
        NewProBnd(i) = ZERO
        NewProAng(i) = ZERO
        NewProDih(i) = ZERO
        NewProUBy(i) = ZERO
        NewProNBC(i) = ZERO
        NewProNBLJ(i) = ZERO
        OldProBnd(i) = ZERO
        OldProAng(i) = ZERO
        OldProDih(i) = ZERO
        OldProUBy(i) = ZERO
        OldProNBC(i) = ZERO
        OldProNBLJ(i) = ZERO
        AvgProBnd(i) = ZERO
        AvgProAng(i) = ZERO
        AvgProDih(i) = ZERO
        AvgProUBy(i) = ZERO
        AvgProNBC(i) = ZERO
        AvgProNBLJ(i) = ZERO
        AvgProBnd2(i) = ZERO
        AvgProAng2(i) = ZERO
        AvgProDih2(i) = ZERO
        AvgProUBy2(i) = ZERO
        AvgProNBC2(i) = ZERO
        AvgProNBLJ2(i) = ZERO
        NewProBndF(i) = ZERO
        NewProAngF(i) = ZERO
        NewProDihF(i) = ZERO
        NewProUByF(i) = ZERO
        NewProNBCF(i) = ZERO
        NewProNBLJF(i) = ZERO
        OldProBndF(i) = ZERO
        OldProAngF(i) = ZERO
        OldProDihF(i) = ZERO
        OldProUByF(i) = ZERO
        OldProNBCF(i) = ZERO
        OldProNBLJF(i) = ZERO
        AvgProBndF(i) = ZERO
        AvgProAngF(i) = ZERO
        AvgProDihF(i) = ZERO
        AvgProUByF(i) = ZERO
        AvgProNBCF(i) = ZERO
        AvgProNBLJF(i) = ZERO
        AvgProBndF2(i) = ZERO
        AvgProAngF2(i) = ZERO
        AvgProDihF2(i) = ZERO
        AvgProUByF2(i) = ZERO
        AvgProNBCF2(i) = ZERO
        AvgProNBLJF2(i) = ZERO
        NewProBndB(i) = ZERO
        NewProAngB(i) = ZERO
        NewProDihB(i) = ZERO
        NewProUByB(i) = ZERO
        NewProNBCB(i) = ZERO
        NewProNBLJB(i) = ZERO
        OldProBndB(i) = ZERO
        OldProAngB(i) = ZERO
        OldProDihB(i) = ZERO
        OldProUByB(i) = ZERO
        OldProNBCB(i) = ZERO
        OldProNBLJB(i) = ZERO
        AvgProBndB(i) = ZERO
        AvgProAngB(i) = ZERO
        AvgProDihB(i) = ZERO
        AvgProUByB(i) = ZERO
        AvgProNBCB(i) = ZERO
        AvgProNBLJB(i) = ZERO
        AvgProBndB2(i) = ZERO
        AvgProAngB2(i) = ZERO
        AvgProDihB2(i) = ZERO
        AvgProUByB2(i) = ZERO
        AvgProNBCB2(i) = ZERO
        AvgProNBLJB2(i) = ZERO
        do j=1,MAXPROTEINS
          NewProProC(i,j) = ZERO
          NewProProLJ(i,j) = ZERO
          NewProProDist(i,j) = ZERO
          OldProProC(i,j) = ZERO
          OldProProLJ(i,j) = ZERO
          OldProProDist(i,j) = ZERO
          AvgProProC(i,j) = ZERO
          AvgProProC2(i,j) = ZERO
          AvgProProLJ(i,j) = ZERO
          AvgProProLJ2(i,j) = ZERO
          NewProProCF(i,j) = ZERO
          NewProProLJF(i,j) = ZERO
          OldProProCF(i,j) = ZERO
          OldProProLJF(i,j) = ZERO
          AvgProProCF(i,j) = ZERO
          AvgProProCF2(i,j) = ZERO
          AvgProProLJF(i,j) = ZERO
          AvgProProLJF2(i,j) = ZERO
          NewProProCB(i,j) = ZERO
          NewProProLJB(i,j) = ZERO
          OldProProCB(i,j) = ZERO
          OldProProLJB(i,j) = ZERO
          AvgProProCB(i,j) = ZERO
          AvgProProCB2(i,j) = ZERO
          AvgProProLJB(i,j) = ZERO
          AvgProProLJB2(i,j) = ZERO
        enddo
        NewProSvnC(i) = ZERO
        NewProSvnLJ(i) = ZERO
        NewProSvnDist(i) = ZERO
        OldProSvnC(i) = ZERO
        OldProSvnLJ(i) = ZERO
        OldProSvnDist(i) = ZERO
        AvgProSvnC(i) = ZERO
        AvgProSvnC2(i) = ZERO
        AvgProSvnLJ(i) = ZERO
        AvgProSvnLJ2(i) = ZERO
        NewProSvnCF(i) = ZERO
        NewProSvnLJF(i) = ZERO
        OldProSvnCF(i) = ZERO
        OldProSvnLJF(i) = ZERO
        AvgProSvnCF(i) = ZERO
        AvgProSvnCF2(i) = ZERO
        AvgProSvnLJF(i) = ZERO
        AvgProGCSLJF2(i) = ZERO
        AvgProGCSLJF(i) = ZERO
        AvgProSvnLJF2(i) = ZERO
        NewProSvnCB(i) = ZERO
        NewProSvnLJB(i) = ZERO
        OldProSvnCB(i) = ZERO
        OldProSvnLJB(i) = ZERO
        AvgProSvnCB(i) = ZERO
        AvgProSvnCB2(i) = ZERO
        AvgProSvnLJB(i) = ZERO
        AvgProSvnLJB2(i) = ZERO
      enddo

      do i=1,MAXSOLUTES
        NewSolBnd(i) = ZERO
        NewSolAng(i) = ZERO
        NewSolDih(i) = ZERO
        NewSolUBy(i) = ZERO
        NewSolNBC(i) = ZERO
        NewSolNBLJ(i) = ZERO
        OldSolBnd(i) = ZERO
        OldSolAng(i) = ZERO
        OldSolDih(i) = ZERO
        OldSolUBy(i) = ZERO
        OldSolNBC(i) = ZERO
        OldSolNBLJ(i) = ZERO
        NewSolBndF(i) = ZERO
        NewSolAngF(i) = ZERO
        NewSolDihF(i) = ZERO
        NewSolUByF(i) = ZERO
        NewSolNBCF(i) = ZERO
        NewSolNBLJF(i) = ZERO
        OldSolBndF(i) = ZERO
        OldSolAngF(i) = ZERO
        OldSolDihF(i) = ZERO
        OldSolUByF(i) = ZERO
        OldSolNBCF(i) = ZERO
        OldSolNBLJF(i) = ZERO
        NewSolBndB(i) = ZERO
        NewSolAngB(i) = ZERO
        NewSolDihB(i) = ZERO
        NewSolUByB(i) = ZERO
        NewSolNBCB(i) = ZERO
        NewSolNBLJB(i) = ZERO
        OldSolBndB(i) = ZERO
        OldSolAngB(i) = ZERO
        OldSolDihB(i) = ZERO
        OldSolUByB(i) = ZERO
        OldSolNBCB(i) = ZERO
        OldSolNBLJB(i) = ZERO
        DoneSolNrg(i) = .false.
        IgnoreSolute(i) = .false.
        NSolVariable(i) = 0
        NSolForwards(i) = 0
        NSolBackwards(i) = 0   
        DualTopologySolute(i) = 0     
        DualTopologySyncTrans(i) = .false.                  
        DualTopologySyncRot(i) = .false.                  
        DualTopologySyncInt(i) = .false.                  
        AvgSolBnd(i) = ZERO           
        AvgSolAng(i) = ZERO           
        AvgSolDih(i) = ZERO           
        AvgSolUBy(i) = ZERO
        AvgSolNBC(i) = ZERO           
        AvgSolNBLJ(i) = ZERO           
        AvgSolBnd2(i) = ZERO           
        AvgSolAng2(i) = ZERO           
        AvgSolDih2(i) = ZERO           
        AvgSolNBC2(i) = ZERO           
        AvgSolNBLJ2(i) = ZERO    
        AvgSolBndF(i) = ZERO           
        AvgSolAngF(i) = ZERO           
        AvgSolDihF(i) = ZERO           
        AvgSolUByF(i) = ZERO
        AvgSolNBCF(i) = ZERO           
        AvgSolNBLJF(i) = ZERO           
        AvgSolBndF2(i) = ZERO           
        AvgSolAngF2(i) = ZERO           
        AvgSolDihF2(i) = ZERO           
        AvgSolNBCF2(i) = ZERO           
        AvgSolNBLJF2(i) = ZERO    
        AvgSolBndB(i) = ZERO           
        AvgSolAngB(i) = ZERO           
        AvgSolDihB(i) = ZERO           
        AvgSolUByB(i) = ZERO
        AvgSolNBCB(i) = ZERO           
        AvgSolNBLJB(i) = ZERO           
        AvgSolBndB2(i) = ZERO           
        AvgSolAngB2(i) = ZERO           
        AvgSolDihB2(i) = ZERO           
        AvgSolNBCB2(i) = ZERO           
        AvgSolNBLJB2(i) = ZERO    
        SolADel(i) = ZERO
        SolRDel(i) = ZERO
        SoluteSoftCore(i) = .false.
        do j=1,MAXSOLUTEATOMS
           SoluteAtomSoftCore(i,j) = .false.
        enddo
        do j=1,MAXSOLUTES
          NewSolSolC(i,j) = ZERO
          NewSolSolLJ(i,j) = ZERO
          NewSolSolDist(i,j) = ZERO
          OldSolSolC(i,j) = ZERO
          OldSolSolLJ(i,j) = ZERO
          OldSolSolDist(i,j) = ZERO
          AvgSolSolC(i,j) = ZERO
          AvgSolSolLJ(i,j) = ZERO
          AvgSolSolC2(i,j) = ZERO
          AvgSolSolLJ2(i,j) = ZERO
          NewSolSolCF(i,j) = ZERO
          NewSolSolLJF(i,j) = ZERO
          OldSolSolCF(i,j) = ZERO
          OldSolSolLJF(i,j) = ZERO
          AvgSolSolCF(i,j) = ZERO
          AvgSolSolLJF(i,j) = ZERO
          AvgSolSolCF2(i,j) = ZERO
          AvgSolSolLJF2(i,j) = ZERO
          NewSolSolCB(i,j) = ZERO
          NewSolSolLJB(i,j) = ZERO
          OldSolSolCB(i,j) = ZERO
          OldSolSolLJB(i,j) = ZERO
          AvgSolSolCB(i,j) = ZERO
          AvgSolSolLJB(i,j) = ZERO
          AvgSolSolCB2(i,j) = ZERO
          AvgSolSolLJB2(i,j) = ZERO
        enddo
        do j=1,MAXPROTEINS
          NewProSolC(j,i) = ZERO
          NewProSolLJ(j,i) = ZERO
          NewProSolDist(j,i) = ZERO
          OldProSolC(j,i) = ZERO
          OldProSolLJ(j,i) = ZERO
          OldProSolDist(j,i) = ZERO
          AvgProSolC(j,i) = ZERO
          AvgProSolLJ(j,i) = ZERO
          AvgProSolC2(j,i) = ZERO
          AvgProSolLJ2(j,i) = ZERO
          NewProSolCF(j,i) = ZERO
          NewProSolLJF(j,i) = ZERO
          OldProSolCF(j,i) = ZERO
          OldProSolLJF(j,i) = ZERO
          AvgProSolCF(j,i) = ZERO
          AvgProSolLJF(j,i) = ZERO
          AvgProSolCF2(j,i) = ZERO
          AvgProSolLJF2(j,i) = ZERO
          NewProSolCB(j,i) = ZERO
          NewProSolLJB(j,i) = ZERO
          OldProSolCB(j,i) = ZERO
          OldProSolLJB(j,i) = ZERO
          AvgProSolCB(j,i) = ZERO
          AvgProSolLJB(j,i) = ZERO
          AvgProSolCB2(j,i) = ZERO
          AvgProSolLJB2(j,i) = ZERO
        enddo
        do j=1,MAXSOLVARIABLES
          SolVariableTyp(i,j) = 0
          SolVariableAtm(i,j) = 0
          SolVariableVal(i,j,1) = ZERO
          SolVariableVal(i,j,1) = ZERO
        enddo
        
        NewSolSvnC(i) = ZERO
        NewSolSvnLJ(i) = ZERO
        NewSolSvnDist(i) = ZERO
        OldSolSvnC(i) = ZERO
        OldSolSvnLJ(i) = ZERO
        OldSolSvnDist(i) = ZERO
        AvgSolSvnC(i) = ZERO
        AvgSolSvnC2(i) = ZERO
        AvgSolSvnLJ(i) = ZERO
        AvgSolSvnLJ2(i) = ZERO
        NewSolSvnCF(i) = ZERO
        NewSolSvnLJF(i) = ZERO
        OldSolSvnCF(i) = ZERO
        OldSolSvnLJF(i) = ZERO
        AvgSolSvnCF(i) = ZERO
        AvgSolSvnCF2(i) = ZERO
        AvgSolSvnLJF(i) = ZERO
        AvgSolSvnLJF2(i) = ZERO
        NewSolSvnCB(i) = ZERO
        NewSolSvnLJB(i) = ZERO
        OldSolSvnCB(i) = ZERO
        OldSolSvnLJB(i) = ZERO
        AvgSolSvnCB(i) = ZERO
        AvgSolSvnCB2(i) = ZERO
        AvgSolSvnLJB(i) = ZERO
        AvgSolSvnLJB2(i) = ZERO
      enddo
      
      do i=1,MAXSOLVENTS
        DoneSvnNrg(i) = .false. 
        SvnProb(i) = ONE
        SvnNewProb(i) = ONE
        SvnPrefSol = 0
        NSvnTry(i) = 0
        NSvnAcp(i) = 0
        NDelTry(i) = 0
        NInsTry(i) = 0
        NDelAcp(i) = 0
        NInsAcp(i) = 0
        NTheAcp(i) = 0
        NTheTry(i) = 0
        NGCSTry(i) = 0
        NGCSAcp(i) = 0
        NSamTry(i) = 0
        NSamAcp(i) = 0
        DoneGCSNrg(i) = .false.        
        GCsoluteSoftCore(i) = .false.
      enddo
      SvnProbParam = ZERO
      SvnSumProb = dble(MAXSOLVENTS)
      SvnNewSumProb = SvnSumProb
      SvnMaxProb = ONE
      SvnNewMaxProb = ONE
c     initialise the residue stacks
      newStackDepth = 0
      oldStackDepth = 0
      do i=1,MAXRESSTACKDEPTH
        do j=1,MAXPROTEINS
          oldStackResProC(i,j) = ZERO
          oldStackResProL(i,j) = ZERO
          oldStackResProDist(i,j) = ZERO
          newStackResProC(i,j) = ZERO
          newStackResProL(i,j) = ZERO
          newStackResProDist(i,j) = ZERO
          oldStackResProCF(i,j) = ZERO
          oldStackResProLF(i,j) = ZERO
          newStackResProCF(i,j) = ZERO
          newStackResProLF(i,j) = ZERO
          oldStackResProCB(i,j) = ZERO
          oldStackResProLB(i,j) = ZERO
          newStackResProCB(i,j) = ZERO
          newStackResProLB(i,j) = ZERO
        enddo
        do j=1,MAXSOLUTES
          oldStackResSolC(i,j) = ZERO
          oldStackResSolL(i,j) = ZERO
          oldStackResSolDist(i,j) = ZERO
          newStackResSolC(i,j) = ZERO
          newStackResSolL(i,j) = ZERO
          newStackResSolDist(i,j) = ZERO
          oldStackResSolCF(i,j) = ZERO
          oldStackResSolLF(i,j) = ZERO
          newStackResSolCF(i,j) = ZERO
          newStackResSolLF(i,j) = ZERO
          oldStackResSolCB(i,j) = ZERO
          oldStackResSolLB(i,j) = ZERO
          newStackResSolCB(i,j) = ZERO
          newStackResSolLB(i,j) = ZERO
        enddo
        oldStackResSvnC(i) = ZERO
        oldStackResSvnL(i) = ZERO
        oldStackResSvnDist(i) = ZERO
        newStackResSvnC(i) = ZERO
        newStackResSvnL(i) = ZERO
        newStackResSvnDist(i) = ZERO
        oldStackResSvnCF(i) = ZERO
        oldStackResSvnLF(i) = ZERO
        newStackResSvnCF(i) = ZERO
        newStackResSvnLF(i) = ZERO
        oldStackResSvnCB(i) = ZERO
        oldStackResSvnLB(i) = ZERO
        newStackResSvnCB(i) = ZERO
        newStackResSvnLB(i) = ZERO
      
        oldStackResBnd(i) = ZERO
        oldStackResAng(i) = ZERO
        oldStackResDih(i) = ZERO
        oldStackResUBy(i) = ZERO
        oldStackResNBC(i) = ZERO
        oldStackResNBL(i) = ZERO
        newStackResBnd(i) = ZERO
        newStackResAng(i) = ZERO
        newStackResDih(i) = ZERO
        newStackResUBy(i) = ZERO
        newStackResNBC(i) = ZERO
        newStackResNBL(i) = ZERO
        oldStackResBndF(i) = ZERO
        oldStackResAngF(i) = ZERO
        oldStackResDihF(i) = ZERO
        oldStackResUByF(i) = ZERO
        oldStackResNBCF(i) = ZERO
        oldStackResNBLF(i) = ZERO
        newStackResBndF(i) = ZERO
        newStackResAngF(i) = ZERO
        newStackResDihF(i) = ZERO
        newStackResUByF(i) = ZERO
        newStackResNBCF(i) = ZERO
        newStackResNBLF(i) = ZERO
        oldStackResBndB(i) = ZERO
        oldStackResAngB(i) = ZERO
        oldStackResDihB(i) = ZERO
        oldStackResUByB(i) = ZERO
        oldStackResNBCB(i) = ZERO
        oldStackResNBLB(i) = ZERO
        newStackResBndB(i) = ZERO
        newStackResAngB(i) = ZERO
        newStackResDihB(i) = ZERO
        newStackResUByB(i) = ZERO
        newStackResNBCB(i) = ZERO
        newStackResNBLB(i) = ZERO

        oldStackResNrg(i) = ZERO
        newStackResNrg(i) = ZERO
        oldStackResNrgF(i) = ZERO
        newStackResNrgF(i) = ZERO
        oldStackResNrgB(i) = ZERO
        newStackResNrgB(i) = ZERO
      enddo

      coordDepth = 0
      solCoordDepth = 0
c      GCSCavity = 0
      NSolGroups = 0
      do i=1,MAXSOLUTES
        NGroupSols(i) = 0
      enddo

      probProtein = 0
      probSolute = 0
      probSolvent = 0
      probTitrate = 0
      probTheta = 0
      probVolume = 0
      probInsertion = 0
      probDeletion = 0
      probGCSolute = 0
      pProtein = 0
      pVolume = 0
      pSolute = 0
      pSolvent = 0
      pTitrate = 0
      pInsertion = 0
      pTheta = 0
      pDeletion = 0
      pGCSolute = 0
      printMoveEvery = 0
                       
      AVGTOTALENERGY = ZERO
      AVGTOTALENERGY2 = ZERO
      AVGTOTALENERGYF = ZERO
      AVGTOTALENERGYF2 = ZERO
      AVGTOTALENERGYB = ZERO
      AVGTOTALENERGYB2 = ZERO
      AVGBACKENERGY = ZERO  
      AVGBACKENERGY2 = ZERO
      AVGFORWARDSENERGY = ZERO
      AVGFORWARDSENERGY = ZERO
      NSTEPS = 0
      AvgVolume = ZERO
      AvgVolume2 = ZERO

      AvgGCSCavity = ZERO

c     initialise everything to do with the extra energies
      call initialiseExtraEnergies(.false.)
            
      gotInitialTotalEnergy = .false.

c     Initialise Labels
      do i=1,MAXID
        LabelID(i) = 0
        LabelType(i) = 0
        LabelMolecule(i) = 0
        LabelResidueName(i) = NULLNAME
        LabelResidueNumber(i) = 0
        LabelName(i) = NULLNAME
      enddo
c     Soft Core Parameters
c     SG, change default values
      CoulSoftPower = 1
      LJSoftPower = 1
      GBSoftPower = 1
      SoftPower = 6.0d0
      SoftPowerInv = ONE / SoftPower
      SoftPowerInv6 = 6.0d0 / SoftPower
      SoftPower2 = SoftPower / 2.0
      DeltaCore = 0.5
      DeltaCoreCoulL = 2.0
      DeltaCoreCoul = DeltaCoreCoulL**SoftPower
      SoftType = SOFTOLD

c     JM initialise surface area module            
      doSA = .false.            
      probe = 0            
      nsurfpoints = 0            
            
      do i=1,MAXSATYPES            
        GotSurface(i,1) = .false.            
        GotSurface(i,2) = .false.            
        SAType(i) = NULLCLJ            
        SARadius(i) = ZERO            
        SASurftens(i) = ZERO            
      enddo            
            
      OldTotalSurface = ZERO            
      OldTotalSurfaceF = ZERO            
      OldTotalSurfaceB = ZERO            
      OldSAEnergy = ZERO            
      OldSAEnergyF = ZERO            
      OldSAEnergyB = ZERO            
            
c         This takes a while to initialise !!            
      do i=1,MAXSOLUTES            
        SoluteredoSA(i) = .false.            
c        OldSoluteSurf(i) = ZERO            
c        NewSoluteSurf(i) = ZERO            
c        do j=1,MAXSOLUTEATOMS            
c          OldSolAtomSurf(i,j) = ZERO            
c          NewSolAtomSurf(i,j) = ZERO            
c        enddo            
      enddo            
      do i=1,MAXPROTEINS            
        do j=1,MAXRESIDUES            
          ResidueredoSA(i,j) = .false.            
c          OldResidueSurf(i,j) = ZERO            
c          NewResidueSurf(i,j) = ZERO            
c          do k=1,MAXATOMS            
c            OldProtAtomSurf(i,k) = ZERO            
c            NewProtAtomSurf(i,k) = ZERO            
c          enddo            
        enddo            
      enddo            
            
      do i=1,MAXSOLVENTS            
        SolventredoSA(i) = .false.            
      enddo            
            
c     JM initialise gb module            
      doGB = .false.            
      do i=1,MAXGBTYPES            
        GotGB(i,1) = .false.            
        GotGB(i,2) = .false.            
        GBType(i) = NULLCLJ           
        GBIBorn(i) = ZERO            
        GBScale(i) = ZERO            
      enddo            
            
      countGBsolute = 0

      end                         
