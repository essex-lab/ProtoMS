
c     Number of simulation iterations
      integer niter
c     Number of RETI blocks
      integer nretiblock
c     Number of PT blocks
      integer nptblock
c     RETI slave python script
      character*255 slave
c     output directory
      character*255 outdir
c     temperature
      double precision retitemp
c     Number of lambda values
      integer nlambda
c     Number of temperatures for each lambda
      integer ntemp(MAXLAMBDAS)
c     Number of requires slaves
      integer nminslaves
c     Lambda values
      double precision lambdaval(MAXLAMBDAS)      
c     Lambda values in the simulation
      double precision lambdacurr(MAXLAMBDAS)
c     Temperatures
      double precision pttemperature(MAXLAMBDAS,MAXTEMPS)
c     Temperatures in the simulation
      double precision pttempcurr(MAXLAMBDAS,MAXTEMPS)
c     Beta factor
      double precision mt_beta_t 
c     Restart filenames
      character*255 rstfiles(MAXLAMBDAS,MAXTEMPS)
c     RETI energy filenames
      character*255 retifiles(MAXLAMBDAS)
c     PT energy filenames
      character*255 energyfiles(MAXLAMBDAS,MAXTEMPS)
c     PT trajectory files
      character*255 pttrajfiles(MAXLAMBDAS)
c     PT file pointer
      integer pttrajstream(MAXLAMBDAS)
c     PT energy dist files
      character*255 enedistfiles(MAXLAMBDAS)
c     PT energy dist pointer
      integer enediststream(MAXLAMBDAS)
c     Attempted number of RETI swaps
      integer nretitry
c     Successful RETI swaps
      integer nretiswap
c     Attempted number of PT swaps
      integer npttry(MAXLAMBDAS)
c     Successful PT swaps
      integer nptswap(MAXLAMBDAS)
      
      common /RETIVAR/ lambdaval,lambdacurr,retitemp,mt_beta_t,niter,nretiblock,nlambda,nretitry,nretiswap,slave,outdir,
     .                   rstfiles,retifiles
      common /RETIVAR2/ pttemperature,pttempcurr,nptblock,ntemp,nminslaves,npttry,nptswap,pttrajstream,
     .                    enediststream,energyfiles,pttrajfiles,enedistfiles

      save /RETIVAR/
      save /RETIVAR2/
