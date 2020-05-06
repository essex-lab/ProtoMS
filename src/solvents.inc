
c   This file contains all of the variables that hold the solvent coordinates and parameters

c     the number of solvents loaded...
      integer NSolvents
      integer NGCsolutes
c     the number of solvent atoms
      integer NSvnAtoms(MAXSOLVENTS)
      integer NGCSAtoms(1000)
c     the cartesian coordinates (the solvent is rigid)
      double precision ZSvnCoords(MAXSOLVENTS,MAXSOLVENTATOMS,3)
      double precision ZGCSCoords(2000,MAXSOLUTEATOMS,3)
c     the mean centre, major minor and perpendicular axis of each solvent
      double precision SvnCentre(MAXSOLVENTS,3)
      double precision GCSCentre(MAXSOLVENTS,3)
c     the old solvent centre - used for volume moves
      double precision OldSvnCentre(MAXSOLVENTS,3)
c      double precision SvnMajAxis(MAXSOLVENTS,3)
c      double precision SvnMinAxis(MAXSOLVENTS,3)
c      double precision SvnPerAxis(MAXSOLVENTS,3)
c     the solute from which preferential sampling is performed
      integer SvnPrefSol
c     variables needed to do preferential sampling...
c     The current probability of moving a particular solvent molecule
      double precision SvnProb(MAXSOLVENTS)
c     The new probability of moving a particular solvent molecule
      double precision SvnNewProb(MAXSOLVENTS)
c     The value of the maximum probability of moving a solvent molecule
      double precision SvnMaxProb,SvnNewMaxProb
c     The sum of the solvent probabilities
      double precision SvnSumProb,SvnNewSumProb
c     The parameter that controls the preferential sampling - WKC
      double precision SvnProbParam
c     how much to translate each solvent by (A)
      double precision SvnADel(MAXSOLVENTS)
      double precision GCSADel(MAXSOLVENTS)
c     how much to rotate each solvent by (rads)
      double precision SvnRDel(MAXSOLVENTS)
      double precision GCSRDel(MAXSOLVENTS)
c     the parameters
      integer ZSvnPar(MAXSOLVENTS,MAXSOLVENTATOMS)
      integer ZGCSPar(MAXSOLVENTS,50)
c     the name of each solvent
      character*4 SvnName(MAXSOLVENTS)
      character*4 GCSname(MAXSOLVENTS)
c     the residue number of the solvent
      integer SvnRes(MAXSOLVENTS)
      integer GCSRes(MAXSOLVENTS)
c     the name of each atom
      character*4 SvnAtomName(MAXSOLVENTS,MAXSOLVENTATOMS)
      character*4 GCSAtomName(MAXSOLVENTS,MAXSOLUTEATOMS)
c     whether or not we have calculated the energy of each solvent
      logical DoneSvnNrg(MAXSOLVENTS)      
c     the theta parameter for each gcs
      logical DoneGCSNrg(MAXSOLVENTS)
      logical GCsoluteSoftCore(MAXSOLVENTS)
      double precision thetaparam(MAXSOLVENTS)
c     The optimization flag for each solvent
      integer SvnOptimization(MAXSOLVENTS)


      common /DSVNDATA/ ZSvnCoords,SvnCentre,OldSvnCentre,
     .                  SvnADel,SvnRDel,SvnProb,SvnNewProb,
     .                  SvnMaxProb,SvnNewMaxProb,SvnSumProb,SvnNewSumProb,
     .                  SvnProbParam,ZGCSCoords,GCSADel,GCSRDel,GCSCentre,SvnOptimization
      common /ISVNDATA/ NSvnAtoms,ZSvnPar,ZGCSPar,NSolvents,SvnRes,SvnPrefSol
      common /CSVNDATA/ SvnName,SvnAtomName,GCSAtomName,GCSname
      common /LSVNDATA/ DoneSvnNrg,DoneGCSNrg
      common /NSVNDATA/ NGCSAtoms
      common /GSVNDATA/ thetaparam,NGCsolutes,GCSRes,GCsoluteSoftCore

      save /DSVNDATA/
      save /ISVNDATA/
      save /CSVNDATA/
      save /LSVNDATA/
      save /NSVNDATA/
      save /GSVNDATA/
