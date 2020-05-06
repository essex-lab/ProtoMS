c
c     File created by Julien Michel November 2004
c     Include files for Surface Area calculations with ProtoMS
c
c     The minimum size of a radius for Hydrogens
      double precision MINRADH
      parameter(MINRADH=1.10D+00)
c     Whether or not to do Surface Area calculations
      logical doSA
c     Whether or not this atom contributes to the non polar energy
      logical SolAtSAcount(MAXSOLUTES,MAXSOLUTEATOMS)
      logical ProtAtSAcount(MAXPROTEINS,MAXATOMS)
      logical SvnAtSAcount(MAXSOLVENTS,MAXSOLVENTATOMS)
c     Whether or not the surface area of this atom may have changed
      logical SolAtredoSA(MAXSOLUTES,MAXSOLUTEATOMS)
      logical ProtAtredoSA(MAXPROTEINS,MAXATOMS)
      logical SvnAtredoSA(MAXSOLVENTS,MAXSOLVENTATOMS)
c    Whether or not this residue needs to be done after a move
      logical SoluteredoSA(MAXSOLUTES)
      logical ResidueredoSA(MAXPROTEINS,MAXRESIDUES)
      logical SolventredoSA(MAXSOLVENTS)

c     If this residue has just been moved
      integer SAMovedResidue(MAXPROTEINS,MAXRESIDUES)
c     If this solute has just been moved
      integer SAMovedSolute(MAXSOLUTES)
c     If this solvent has just been moved
      integer SAMovedSolvent(MAXSOLVENTS)
      common /ISURFAREA/ SAMovedResidue,SAMovedSolute,SAMovedSolvent

      save /ISURFAREA/
c
c     Holds the radii used for surface area calculations
c
      integer MAXSATYPES
      parameter(MAXSATYPES=200)
      character*2 SAType(MAXSATYPES)
      logical GotSurface(MAXSATYPES,2)
      double precision SARadius(MAXSATYPES)
      double precision SASurftens(MAXSATYPES)
      double precision SurfRadDB(MAXCLJ)
      double precision SurfTensionDB(MAXCLJ)
c    This is for perturbed atoms
      double precision VariableSurfRadius(MAXVARIABLEPARAM,5)
      double precision VariableSurfTension(MAXVARIABLEPARAM,5)

c      double precision SolAtRadius(MAXSOLUTES,MAXSOLUTEATOMS)
c      double precision ProtAtRadius(MAXPROTEINS,MAXATOMS)
c      double precision SvnAtRadius(MAXSOLVENTS,MAXSOLVENTATOMS)
c
c
c     Holds the Surface Area of every atom
c
      double precision OldTotalSurface
      double precision NewTotalSurface
      double precision OldTotalSurfaceF
      double precision NewTotalSurfaceF
      double precision OldTotalSurfaceB
      double precision NewTotalSurfaceB
      double precision AvgTotalSurface
      double precision AvgTotalSurfaceF
      double precision AvgTotalSurfaceB
      double precision AvgTotalSurface2
      double precision AvgTotalSurfaceF2
      double precision AvgTotalSurfaceB2
      double precision OldProtAtomSurf(MAXPROTEINS,MAXATOMS)
      double precision NewProtAtomSurf(MAXPROTEINS,MAXATOMS)
      double precision OldProtAtomSurfF(MAXPROTEINS,MAXATOMS)
      double precision NewProtAtomSurfF(MAXPROTEINS,MAXATOMS)
      double precision OldProtAtomSurfB(MAXPROTEINS,MAXATOMS)
      double precision NewProtAtomSurfB(MAXPROTEINS,MAXATOMS)
      double precision OldResidueSurf(MAXPROTEINS,MAXRESIDUES)
      double precision NewResidueSurf(MAXPROTEINS,MAXRESIDUES)
      double precision OldResidueSurfF(MAXPROTEINS,MAXRESIDUES)
      double precision NewResidueSurfF(MAXPROTEINS,MAXRESIDUES)
      double precision OldResidueSurfB(MAXPROTEINS,MAXRESIDUES)
      double precision NewResidueSurfB(MAXPROTEINS,MAXRESIDUES)
      double precision OldProteinSurf(MAXPROTEINS)
      double precision NewProteinSurf(MAXPROTEINS)
      double precision OldProteinSurfF(MAXPROTEINS)
      double precision NewProteinSurfF(MAXPROTEINS)
      double precision OldProteinSurfB(MAXPROTEINS)
      double precision NewProteinSurfB(MAXPROTEINS)
      double precision OldSolAtomSurf(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision NewSolAtomSurf(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision OldSolAtomSurfF(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision NewSolAtomSurfF(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision OldSolAtomSurfB(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision NewSolAtomSurfB(MAXSOLUTES,MAXSOLUTEATOMS)
      double precision OldSoluteSurf(MAXSOLUTES)
      double precision NewSoluteSurf(MAXSOLUTES)
      double precision OldSoluteSurfF(MAXSOLUTES)
      double precision NewSoluteSurfF(MAXSOLUTES)
      double precision OldSoluteSurfB(MAXSOLUTES)
      double precision NewSoluteSurfB(MAXSOLUTES)
      double precision OldSvnAtomSurf(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision NewSvnAtomSurf(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision OldSvnAtomSurfF(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision NewSvnAtomSurfF(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision OldSvnAtomSurfB(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision NewSvnAtomSurfB(MAXSOLVENTS,MAXSOLVENTATOMS)
      double precision OldSolventSurF(MAXSOLVENTS)
      double precision NewSolventSurf(MAXSOLVENTS)
      double precision OldSolventSurFF(MAXSOLVENTS)
      double precision NewSolventSurfF(MAXSOLVENTS)
      double precision OldSolventSurFB(MAXSOLVENTS)
      double precision NewSolventSurfB(MAXSOLVENTS)
c     Cutoff for speed up
      double precision SurfCut,SurfCut2
c     Energy components
      double precision OldSAEnergy
      double precision OldSAEnergyf
      double precision OldSAEnergyb
      double precision NewSAEnergy
      double precision NewSAEnergyf
      double precision NewSAEnergyb
      double precision AvgSAEnergy
      double precision AvgSAEnergyF
      double precision AvgSAEnergyB
      double precision AvgSAEnergy2
      double precision AvgSAEnergyF2
      double precision AvgSAEnergyB2

      common /LSURFAREA/ doSA,SolAtSAcount,ProtAtSAcount,SvnAtSAcount,
     .                   SolAtredoSA,ProtAtredoSA,SvnAtredoSA,SoluteredoSA,ResidueredoSA,
     .                   SolventredoSA,GotSurface

      common /DSURFAREA/ OldTotalSurface,NewTotalSurface,
     .  OldTotalSurfaceF,OldTotalSurfaceB,
     .  NewTotalSurfaceF,NewTotalSurfaceB,
     .  AvgTotalSurface,AvgTotalSurfaceF,AvgTotalSurfaceB,
     .  AvgTotalSurface2,AvgTotalSurfaceF2,AvgTotalSurfaceB2,
     .  OldProtAtomSurf,NewProtAtomSurf,OldProtAtomSurfF,NewProtAtomSurfF,
     .  OldProtAtomSurfB,NewProtAtomSurfB,
     .  OldSolAtomSurf,NewSolAtomSurf,OldSolAtomSurfF,NewSolAtomSurfF,
     .  OldSolAtomSurfB,NewSolAtomSurfB, 
     .  OldSvnAtomSurf,NewSvnAtomSurf,     
     .  OldSvnAtomSurfF,NewSvnAtomSurfF,
     .  OldSvnAtomSurfB,NewSvnAtomSurfB,
     .  OldProteinSurf,NewProteinSurf,OldProteinSurfF,NewProteinSurfF,
     .  OldProteinSurfB,NewProteinSurfB, 
     .  OldSoluteSurf,NewSoluteSurf,OldSoluteSurfF,NewSoluteSurfF,
     .  OldSoluteSurfB,NewSoluteSurfB,
     .  OldSolventSurf,NewSolventSurf,
     .  OldSolventSurfF,NewSolventSurfF,
     .  OldSolventSurfB,NewSolventSurfB,
     .  OldResidueSurf,NewResidueSurf,OldResidueSurfF,NewResidueSurfF,
     .  OldResidueSurfB,NewResidueSurfB,
     .  SurfCut,SurfCut2,
     .  SARadius,SASurftens,SurfRadDB,SurfTensionDB,VariableSurfRadius,
     .  VariableSurfTension

      common /DSANRG/ OldSAEnergy,OldSAEnergyf,OldSAEnergyb,
     .  NewSAEnergy,NewSAEnergyf,NewSAEnergyb,AvgSAEnergy,
     .  AvgSAEnergyF,AvgSAEnergyB,AvgSAEnergy2,AvgSAEnergyF2,
     .  AvgSAEnergyB2

      common /CSURFAREA/ SAType
      
      save /LSURFAREA/
      save /DSURFAREA/
      save /DSANRG/
      save /CSURFAREA/

C#####################################################
C Include file for shrake and rupley. It list the three 
C data blocks holding points on the surface of a sphere
C centered on the origin with a unit radius.
C#####################################################

c The number of surface elements
c should be 12,42,162 or 642, the higher, the more accurate the calculation
c is determined by SAquality
      integer nsurfpoints
c can be 1 (12 pts), 2 (42 pts), 3 (162pts), 4 (642pts) 
      integer SAquality
      common /ISHRAKE/ nsurfpoints,SAquality

c Constants
      double precision PIPI
      parameter(PIPI=3.14159265358979323846264D+00)
      double precision FOURPI
      parameter(FOURPI=PIPI*4.0D+00)

c The radius of the probe
      double precision probe
      common /DSHRAKE/ probe

c The coordinates of these points
c the coordinates are in shrakeinit.inc
      double precision points12(36)
      double precision points42(126)
      double precision points162(486)
      double precision points642(1926)

      common /shrake_points/ points12, points42, points162, 
     .                       points642

      save /shrake_points/


      logical DualredoSA
     
      common /LDUALSASA/ DualredoSA
      
