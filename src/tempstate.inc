
      integer tGCSCavity
      double precision tBoxDimension(3),tBoxOrigin(3)
      double precision tZCoords(MAXPROTEINS,MAXATOMS,3)
      double precision tSolDummy(MAXSOLUTES,3,3)
      double precision tSolGeomCentre(MAXSOLUTES,3)
      double precision tZSolCoords(MAXSOLUTES,MAXSOLUTEATOMS,3)
      double precision tSvnCentre(MAXSOLVENTS,3)
      double precision tZSvnCoords(MAXSOLVENTS,MAXSOLVENTATOMS,3)
      double precision tGCSCentre(MAXSOLVENTS,3)
      double precision tZGCSCoords(MAXSOLVENTS,MAXSOLVENTATOMS,3)
      double precision tthetaparam(MAXSOLVENTS)
      
      common /tstate/ tBoxDimension,tBoxOrigin,tZCoords,tSolDummy,tSolGeomCentre,tZSolCoords,tSvnCentre,
     .                tZSvnCoords,tGCSCentre,tZGCSCoords,tthetaparam,tGCSCavity