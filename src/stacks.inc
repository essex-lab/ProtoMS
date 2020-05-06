
c   This file contains the variables used in the coordinate and residue stacks
       

c     the residue energy stack
      double precision oldStackResProC(MAXRESSTACKDEPTH,MAXPROTEINS)
      double precision oldStackResProL(MAXRESSTACKDEPTH,MAXPROTEINS)
      double precision oldStackResProDist(MAXRESSTACKDEPTH,MAXPROTEINS)
      double precision newStackResProC(MAXRESSTACKDEPTH,MAXPROTEINS)
      double precision newStackResProL(MAXRESSTACKDEPTH,MAXPROTEINS)
      double precision newStackResProDist(MAXRESSTACKDEPTH,MAXPROTEINS)
      double precision oldStackResSolC(MAXRESSTACKDEPTH,MAXSOLUTES)
      double precision oldStackResSolL(MAXRESSTACKDEPTH,MAXSOLUTES)
      double precision oldStackResSolDist(MAXRESSTACKDEPTH,MAXSOLUTES)
      double precision newStackResSolC(MAXRESSTACKDEPTH,MAXSOLUTES)
      double precision newStackResSolL(MAXRESSTACKDEPTH,MAXSOLUTES)
      double precision newStackResSolDist(MAXRESSTACKDEPTH,MAXSOLUTES)
      double precision oldStackResSvnC(MAXRESSTACKDEPTH)
      double precision oldStackResSvnL(MAXRESSTACKDEPTH)
      double precision oldStackResSvnDist(MAXRESSTACKDEPTH)
      double precision newStackResSvnC(MAXRESSTACKDEPTH)
      double precision newStackResSvnL(MAXRESSTACKDEPTH)
      double precision newStackResSvnDist(MAXRESSTACKDEPTH)
      double precision oldStackResBnd(MAXRESSTACKDEPTH)
      double precision oldStackResAng(MAXRESSTACKDEPTH)
      double precision oldStackResDih(MAXRESSTACKDEPTH)
      double precision oldStackResUBy(MAXRESSTACKDEPTH)
      double precision oldStackResNBC(MAXRESSTACKDEPTH)
      double precision oldStackResNBL(MAXRESSTACKDEPTH)
      double precision newStackResBnd(MAXRESSTACKDEPTH)
      double precision newStackResAng(MAXRESSTACKDEPTH)
      double precision newStackResDih(MAXRESSTACKDEPTH)
      double precision newStackResUBy(MAXRESSTACKDEPTH)
      double precision newStackResNBC(MAXRESSTACKDEPTH)
      double precision newStackResNBL(MAXRESSTACKDEPTH)
c       forwards residue energy stack
      double precision oldStackResProCF(MAXRESSTACKDEPTH,MAXPROTEINS)
      double precision oldStackResProLF(MAXRESSTACKDEPTH,MAXPROTEINS)
      double precision newStackResProCF(MAXRESSTACKDEPTH,MAXPROTEINS)
      double precision newStackResProLF(MAXRESSTACKDEPTH,MAXPROTEINS)
      double precision oldStackResSolCF(MAXRESSTACKDEPTH,MAXSOLUTES)
      double precision oldStackResSolLF(MAXRESSTACKDEPTH,MAXSOLUTES)
      double precision newStackResSolCF(MAXRESSTACKDEPTH,MAXSOLUTES)
      double precision newStackResSolLF(MAXRESSTACKDEPTH,MAXSOLUTES)
      double precision oldStackResSvnCF(MAXRESSTACKDEPTH)
      double precision oldStackResSvnLF(MAXRESSTACKDEPTH)
      double precision newStackResSvnCF(MAXRESSTACKDEPTH)
      double precision newStackResSvnLF(MAXRESSTACKDEPTH)
      double precision oldStackResBndF(MAXRESSTACKDEPTH)
      double precision oldStackResAngF(MAXRESSTACKDEPTH)
      double precision oldStackResDihF(MAXRESSTACKDEPTH)
      double precision oldStackResUByF(MAXRESSTACKDEPTH)
      double precision oldStackResNBCF(MAXRESSTACKDEPTH)
      double precision oldStackResNBLF(MAXRESSTACKDEPTH)
      double precision newStackResBndF(MAXRESSTACKDEPTH)
      double precision newStackResAngF(MAXRESSTACKDEPTH)
      double precision newStackResDihF(MAXRESSTACKDEPTH)
      double precision newStackResUByF(MAXRESSTACKDEPTH)
      double precision newStackResNBCF(MAXRESSTACKDEPTH)
      double precision newStackResNBLF(MAXRESSTACKDEPTH)
c       backwards energy residue stack
      double precision oldStackResProCB(MAXRESSTACKDEPTH,MAXPROTEINS)
      double precision oldStackResProLB(MAXRESSTACKDEPTH,MAXPROTEINS)
      double precision newStackResProCB(MAXRESSTACKDEPTH,MAXPROTEINS)
      double precision newStackResProLB(MAXRESSTACKDEPTH,MAXPROTEINS)
      double precision oldStackResSolCB(MAXRESSTACKDEPTH,MAXSOLUTES)
      double precision oldStackResSolLB(MAXRESSTACKDEPTH,MAXSOLUTES)
      double precision newStackResSolCB(MAXRESSTACKDEPTH,MAXSOLUTES)
      double precision newStackResSolLB(MAXRESSTACKDEPTH,MAXSOLUTES)
      double precision oldStackResSvnCB(MAXRESSTACKDEPTH)
      double precision oldStackResSvnLB(MAXRESSTACKDEPTH)
      double precision newStackResSvnCB(MAXRESSTACKDEPTH)
      double precision newStackResSvnLB(MAXRESSTACKDEPTH)
      double precision oldStackResBndB(MAXRESSTACKDEPTH)
      double precision oldStackResAngB(MAXRESSTACKDEPTH)
      double precision oldStackResDihB(MAXRESSTACKDEPTH)
      double precision oldStackResUByB(MAXRESSTACKDEPTH)
      double precision oldStackResNBCB(MAXRESSTACKDEPTH)
      double precision oldStackResNBLB(MAXRESSTACKDEPTH)
      double precision newStackResBndB(MAXRESSTACKDEPTH)
      double precision newStackResAngB(MAXRESSTACKDEPTH)
      double precision newStackResDihB(MAXRESSTACKDEPTH)
      double precision newStackResUByB(MAXRESSTACKDEPTH)
      double precision newStackResNBCB(MAXRESSTACKDEPTH)
      double precision newStackResNBLB(MAXRESSTACKDEPTH)
c       the main residue energies (on the stack!)
      double precision newStackResNrg(MAXRESSTACKDEPTH)
      double precision oldStackResNrg(MAXRESSTACKDEPTH)
      double precision newStackResNrgF(MAXRESSTACKDEPTH)
      double precision newStackResNrgB(MAXRESSTACKDEPTH)
      double precision oldStackResNrgF(MAXRESSTACKDEPTH)
      double precision oldStackResNrgB(MAXRESSTACKDEPTH)
      integer newStackDepth,oldStackDepth
            

      common /ENERGYSTACK1/ newstackresnblb,newstackresnbcb,newstackresdihb,
     . newstackresangb,newstackresbndb,oldstackresnblb,oldstackresnbcb,
     . oldStackResNrgF,oldStackResNrgB,newStackResNrgF,newStackResNrgB,
     . oldstackresdihb,oldstackresangb,oldstackresbndb,newstackressvnlb,
     . oldstackresubyb,oldstackresubyf,newstackresubyb,newstackresubyf,
     . newstackressvncb,oldstackressvnlb,oldstackressvncb,newstackressollb,
     . newstackressolcb,oldstackressollb,oldstackressolcb,newstackresprolb,
     . newstackresprocb,oldstackresprolb,oldstackresprocb,newstackresnblf,
     . newstackresnbcf,newstackresdihf,newstackresangf,newstackresbndf,
     . oldstackresnblf,oldstackresnbcf,oldstackresdihf,oldstackresangf,
     . oldstackresbndf,newstackressvnlf,newstackressvncf,oldstackressvnlf,
     . oldstackressvncf,newstackressollf,newstackressolcf,oldstackressollf,
     . oldstackressolcf,newstackresprolf,newstackresprocf,oldstackresprolf,
     . oldstackresprocf
     
      save /ENERGYSTACK1/

      
      common /ENERGYSTACK2/ oldStackResProC,oldStackResProL,oldStackResProDist,
     .                     newStackResProC,newStackResProL,newStackResProDist,
     .                     oldStackResSolC,oldStackResSolL,oldStackResSolDist,
     .                     newStackResSolC,newStackResSolL,newStackResSolDist,
     .                     oldStackResSvnC,oldStackResSvnL,oldStackResSvnDist,
     .                     newStackResSvnC,newStackResSvnL,newStackResSvnDist,
     . oldStackResBnd,oldStackResAng,oldStackResDih,oldStackResNBC,oldStackResNBL,
     . newStackResBnd,newStackResAng,newStackResDih,newStackResNBC,newStackResNBL,
     . newStackResNrg,oldStackResNrg,newStackDepth,oldStackDepth,
     . oldStackResUBy,newStackResUBy
     
      save /ENERGYSTACK2/
     
c     the residue coordinate stack
      double precision stackResBB(MAXRESCOORDDEPTH,4,3)
      double precision stackResSC(MAXRESCOORDDEPTH,MAXSCATOMS,3)      
      double precision stackResZSC(MAXRESCOORDDEPTH,MAXSCATOMS,3)
      double precision stackResCentre(MAXRESCOORDDEPTH,3)
      double precision stackResRadius(MAXRESCOORDDEPTH)
      integer stackResCoord(MAXRESCOORDDEPTH),stackResProt(MAXRESCOORDDEPTH)
      integer coordDepth
      
      common /RESIDUESTACK/ stackResBB,stackResSC,stackResZSC,stackResCentre,
     .                      stackResRadius,stackResCoord,stackResProt,
     .                      coordDepth

      save /RESIDUESTACK/

c     the solute coordinate stack
      double precision stackSol(MAXSOLCOORDDEPTH,MAXSOLUTEATOMS,3)
      double precision stackSolZ(MAXSOLCOORDDEPTH,MAXSOLUTEATOMS,3)
      double precision stackSolBonds(MAXSOLCOORDDEPTH,MAXSOLUTEBONDS)
      double precision stackSolAngles(MAXSOLCOORDDEPTH,MAXSOLUTEANGLES)
      double precision stackSolDihedrals(MAXSOLCOORDDEPTH,MAXSOLUTEDIHEDRALS)
      double precision stackSolUreyBradleys(MAXSOLCOORDDEPTH,MAXSOLUTEUREYBRADLEYS)
      double precision stackSolDummies(MAXSOLCOORDDEPTH,3,3)
      double precision stackSolCentre(MAXSOLCOORDDEPTH,3)
      double precision stackSolGeomCentre(MAXSOLCOORDDEPTH,3)
      double precision stackSolRadius(MAXSOLCOORDDEPTH)
c      double precision stackSolAxis(MAXSOLCOORDDEPTH,3,3)
      integer stackSolCoord(MAXSOLCOORDDEPTH)
      integer solCoordDepth
c      integer namsol(MAXSOLVENTS)
      common /SOLUTESTACK/ stackSol,stackSolZ,stackSolBonds,stackSolAngles,stackSolDihedrals,stackSolUreyBradleys,
     .                     stackSolCentre,stackSolGeomCentre,stackSolDummies,
     .                     stackSolRadius,solCoordDepth,stackSolCoord      
      save /SOLUTESTACK/

c     the solvent coordinate stack
      double precision stackSvn(MAXSVNCOORDDEPTH,MAXSOLVENTATOMS,3)
      double precision stackSvnCentre(MAXSVNCOORDDEPTH,3)
      double precision stackGCS(MAXSVNCOORDDEPTH,MAXSOLUTEATOMS,3)
      double precision stackGCSCentre(MAXSVNCOORDDEPTH,3)
c      double precision stackSvnAxis(MAXSVNCOORDDEPTH,3,3)
      integer stackSvnCoord(MAXSVNCOORDDEPTH)
      integer stackGCSCoord(MAXSVNCOORDDEPTH)
      integer svnCoordDepth
      integer gcsCoordDepth
      
      common /SOLVENTSTACK/ stackSvn,stackGCS,stackSvnCentre,stackGCSCentre,stackSvnCoord,
     .                        svnCoordDepth,stackGCSCoord,gcsCoordDepth
 
      save /SOLVENTSTACK/
 
c     cavity coordinate stacks
c      integer MAXCAVCOORDDEPTH
c      double precision stackCav(2500,2500,1)
c      double precision stackCav(2500,2500,2)
c      double precision ZSolCoordst(1000,1000,3)
c      double precision stackCav(5000,20,3)
c      double precision stackCavZ(MAXBNDPARAM,MAXSOLUTEATOMS,3) 
c      double precision stackCavCentre(MAXBNDPARAM,3) 
c      integer stackCavCoord(MAXBNDPARAM)
      double precision SvnReser(10000,1000)
      double precision SvnCav(10000,1000)
      double precision SvnNameC(10000)
      double precision SvnNameR(10000)
      double precision GCSNameC(10000)
      double precision GCSNameR(10000)
      integer cavCoordDepth
      integer bnzCoordDepth
      integer aceCoordDepth
      integer pyrCoordDepth
      integer filledCoordDepth
      integer SvnCavity
      integer GCSCavity
      integer GCSReservoir
      integer SvnReservoir
      integer SvnNonInter 
      integer GCSNonInter   
      integer highstack
      integer avHighstack
      integer instHighstack
      integer lowstack
      integer avLowstack
      integer instLowstack


c      common /CAVITYSTACK/ ZSolCoordst,stackCav,stackCavZ,stackCavCentre,
c     .                     stackCavCoord,cavCoordDepth,filledCoordDepthisvnc  
c     
       common /CAVITYSTACK/ cavCoordDepth,filledCoordDepth,SvnReser,SvnCav,SvnCavity,SvnReservoir,SvnNameC,
     .                       GCSNameC,GCSNameR,SvnNameR,SvnNonInter,
     .                       GCSCavity,GCSReservoir,GCSNonInter,highstack,lowstack,bnzCoordDepth,aceCoordDepth,pyrCoordDepth,
     .                       avHighstack,avLowstack,instHighstack,instLowstack
       
       save /CAVITYSTACK/

c      save /CAVITYSTACK/  

c     Grid coordinate stacks

      double precision stackGridX(5000000,3)
      double precision stackGridY(5000000,3)
      double precision stackGridZ(5000000,3)
c      double precision stackBNZX(500000,3)
      double precision stackGridXL(5000000,3)
      double precision stackGridYL(5000000,3)
      double precision stackGridZL(5000000,3)
c      double precision stackBNZY(500000,3)
c      double precision stackBNZZ(500000,3)
c      double precision stackPYRX(500000,3)
c      double precision stackPYRY(500000,3)
c      double precision stackPYRZ(500000,3)
c      double precision stackACEX(500000,3)
c      double precision stackACEY(500000,3)
c      double precision stackACEZ(500000,3)
      double precision stackResX(MAXRESERVOIR,3)
      double precision stackResY(MAXRESERVOIR,3)
      double precision stackresZ(MAXRESERVOIR,3)  
      double precision stackTrueX(MAXRESERVOIR,3)
      double precision stackTrueY(MAXRESERVOIR,3)
      double precision stackTrueZ(MAXRESERVOIR,3)
      double precision stackFalseX(MAXRESERVOIR,3)
      double precision stackFalseY(MAXRESERVOIR,3)
      double precision stackFalseZ(MAXRESERVOIR,3)
      double precision pointX(100,100)
      double precision pointY(100,100)
      double precision pointZ(100,100)
      integer reservoir
       integer trues
       integer falses
       integer test

      common /GRID/ stackGridX,stackGridY,stackGridZ,stackTrueX,stackTrueY,stackTrueZ,stackFalseX,stackFalseY,
     .                stackFalseZ,stackResX,stackResY,stackResZ,reservoir,trues,falses,test,pointX,pointY,pointZ
                                                 
      common /G2/ stackGridXL,stackGridYL,stackGridZL
      save /GRID/
      save /G2/

