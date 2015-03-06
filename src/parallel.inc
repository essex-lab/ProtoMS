
c     MPI ID and number of tasks
      integer myid,ntasks 
      integer mylamid,mytemid,mygcid,myjawsid
      integer nlamtasks,ntemtasks,ngctasks,njawstasks
      integer nlamladders,ntemladders,ngcladders,njawsladders
      integer lamwithlad(MAXLAMBDAS),temwithlad(MAXTEMPS),gcwithlad(MAXADAMS),jawswithlad(MAXJAWS2)

c     lambda-replica exchange
      integer nlretry,nlreswap
      integer lrefreq,nlambdas,mylambda
      character*255 reoutfolder,lrefolder
      double precision lambdaval(MAXLAMBDAS)  
      double precision deltalam
      integer task2lambda(MAXLAMBDAS),lambda2task(MAXLAMBDAS)

c     temperature-replica exchange
      integer ntretry,ntrswap
      integer trefreq,ntemperatures,mytemperature
      character*255 trefolder
      double precision temval(MAXTEMPS)
      integer task2temperature(MAXLAMBDAS,MAXTEMPS),temperature2task(MAXLAMBDAS,MAXTEMPS)

c     multi GCMC
      integer nadams,myadam
      character*255 gcmcfolder
      double precision adamparam(MAXLAMBDAS)

c     multi JAWS-2
      integer njaws2,myjbias
      character*255 jaws2folder
      double precision alljbias(MAXLAMBDAS)
             
      common /MPISTUFF/ myid,ntasks
      save /MPISTUFF/
 
      common /LRE/ lambdaval,deltalam,task2lambda,lambda2task,nlretry,nlreswap,lrefreq,nlambdas,
     .               mylambda,lrefolder,reoutfolder
      save /LRE/

      common /TRE/ tempval,task2temperature,temperature2task,ntretry,ntrswap,trefreq,ntemperatures,mytemperature,
     .             trefolder
      save /TRE/

      common /MGCMC/ nadams,myadam,adamparam,gcmcfolder
      save /MGCMC/

      common /MJAWS/ njaws2,myjbias,alljbias,jaws2folder
      save /MJAWS/
