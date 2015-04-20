
c     MPI ID and number of tasks
      integer myid,ntasks 
      integer mylamid,mytemid,mygcid,myjawsid
c      integer nlamtasks,ntemtasks,ngctasks,njawstasks
      integer nlamladders,ntemladders,ngcladders,njawsladders
      double precision lamwithlad(MAXLAMBDAS),temwithlad(MAXTEMPS),gcwithlad(MAXADAMS),jawswithlad(MAXJAWS2)

c     lambda-replica exchange
      integer nlretry,nlreswap,mylamcom
      integer lrefreq,nlambdas,mylambda,mylamladder
      character*255 reoutfolder,lrefolder,myoutfolder
      double precision lambdaval(MAXLAMBDAS)  
      double precision deltalam
      integer task2lambda(MAXLAMBDAS),lambda2task(MAXLAMBDAS)

c     temperature-replica exchange
      integer ntretry,ntreswap,mytemcom
      integer trefreq,ntemperatures,mytemperature,mytemladder
      character*255 trefolder
      double precision temval(MAXTEMPS)
      integer task2temperature(MAXTEMPS),temperature2task(MAXTEMPS)
      double precision pt_beta_t(MAXLAMBDAS,MAXTEMPS)   

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
 
      common /LRE/ lambdaval,lamwithlad,deltalam,task2lambda,lambda2task,nlretry,nlreswap,lrefreq,nlambdas,
     .               mylambda,mylamcom,mylamladder,mylamid,nlamladders,lrefolder,reoutfolder,myoutfolder
      save /LRE/

      common /TRE/ pt_beta_t,temval,temwithlad,task2temperature,temperature2task,ntretry,ntreswap,trefreq,ntemperatures,mytemperature,
     .               mytemcom,mytemladder,mytemid,ntemladders,trefolder
      save /TRE/

      common /MGCMC/ gcwithlad,nadams,myadam,mygcid,ngcladders,adamparam,gcmcfolder
      save /MGCMC/

      common /MJAWS/ jawswithlad,njaws2,myjbias,myjawsid,njawsladders,alljbias,jaws2folder
      save /MJAWS/
