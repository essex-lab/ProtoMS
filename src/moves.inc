
c  This file holds the variables used to control and record the moves

c     the pattern of moves in the simulation
      integer probTotal,probSolute,probSolvent,probProtein,probTitrate
      integer probVolume,probSvnSwap,probRotamer
      integer probDeletion,probInsertion,probGCSolute,probSample,probLambda,probSolTheta
      integer pSolute,pSolvent,pProtein,pTitrate,pVolume,pSvnSwap,pRotamer,pDeletion,pInsertion,pGCSolute,pSample,pLambda
      integer printMoveEvery,probTheta,pTheta,pSolTheta
c     the number of attempted and accepted moves for each moeity
      integer NResTry(MAXPROTEINS,MAXRESIDUES),NResAcp(MAXPROTEINS,MAXRESIDUES)
      integer NSolTry(MAXSOLUTES,MAXSOLUTERESIDUES),NSolAcp(MAXSOLUTES,MAXSOLUTERESIDUES),NInsTry(MAXSOLVENTS),NInsAcp(MAXSOLVENTS)
      integer NSvnTry(MAXSOLVENTS),NSvnAcp(MAXSOLVENTS),NDelTry(MAXSOLVENTS),NDelAcp(MAXSOLVENTS)
      integer NGCSTry(MAXSOLVENTS),NGCSAcp(MAXSOLVENTS),NTheAcp(MAXSOLVENTS),
     .          NTheTry(MAXSOLVENTS),NSamTry(MAXSOLVENTS),NSamAcp(MAXSOLVENTS) 
      integer NSolTheTry(MAXSOLUTES),NSolTheAcp(MAXSOLUTES)     
 
c     Multiple level Ewald moves
      integer NEwaldTry,NEwaldAcp
           
      common /MOVEVALS/ ProbTotal,ProbTitrate,ProbProtein,ProbSolute,ProbDeletion,ProbInsertion,ProbGCSolute,ProbSample,
     .          ProbSolvent,ProbVolume,ProbSvnSwap,ProbRotamer,pProtein,pSolute,pVolume,pDeletion,pInsertion,pGCSolute,
     .          pSample,pLambda,pSolvent,pTitrate,pSvnSwap,pRotamer,printMoveEvery,NTheAcp,NTheTry,NSamTry,NSamAcp,ProbLambda,
     .          NResTry,NResAcp,NSolTry,NSolAcp,NSvnTry,NSvnAcp,NInsTry,NInsAcp,NDelAcp,NDelTry,NGCSTry,NGCSAcp,NEwaldTry,NEwaldAcp,
     .          ProbTheta,pTheta,probSolTheta,pSolTheta,NSolTheTry,NSolTheAcp

      save /MOVEVALS/
      
