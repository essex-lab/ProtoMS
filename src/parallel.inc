
c     MPI ID and number of tasks
      integer myid,ntasks

c     Maximum allowed number of lambdas
      integer MAXLAMBDAS
      parameter (MAXLAMBDAS=50)
c     Maximum allowed number of temperatures
      integer MAXTEMPS
      parameter (MAXTEMPS=50)  

c     lambda-replica exchange
      integer nlretry,nlreswap
      integer lrefreq,nlambdas,mylambda
      character*255 reoutfolder,lrefolder
      double precision lambdaval(MAXLAMBDAS)  
      integer task2lambda(MAXLAMBDAS),lambda2task(MAXLAMBDAS)

c     temperature-replica exchange
      integer ntretry,ntrswap
      integer trefreq,ntemperatures(MAXLAMBDAS),mytemperature
      character*255 trefolder
      double precision tempval(MAXLAMBDAS,MAXTEMPS)
      integer task2temperature(MAXLAMBDAS,MAXTEMPS),temperature2task(MAXLAMBDAS,MAXTEMPS)

c     multi GCMC
      integer nadams,myadam
      character*255 gcmcfolder
      double precision adamparam(MAXLAMBDAS)
             
      common /MPISTUFF/ myid,ntasks
      save /MPISTUFF/
 
      common /LRE/ lambdaval,task2lambda,lambda2task,nlretry,nlreswap,lrefreq,nlambdas,mylambda,lrefolder,reoutfolder
      save /LRE/

      common /TRE/ tempval,task2temperature,temperature2task,ntretry,ntrswap,trefreq,ntemperatures,mytemperature,
     .             trefolder
      save /TRE/

      common /MGCMC/ nadams,myadam,adamparam,gcmcfolder
      save /MGCMC/
