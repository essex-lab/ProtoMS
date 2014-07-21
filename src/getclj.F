      integer function getCLJ(ipar0,ipar1)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'surface.inc'
      include 'gb.inc'
c##############################################################
c
c     This function adds the variable clj parameter
c     that mutates ipar0 at lam=ZERO, to ipar1 at lam=1.0
c
c     It then returns the negative of the index to this
c     variable parameter
c
c     The array has 5 spaces, ordered by;
c
c     1 = Value at lambda = 0
c     2 = Value at lambda = Lambda
c     3 = Value at lambda = 1.0
c     4 = Value at lambda = backwards (lambdaB - for singlecopy zwanzig free energy)
c     5 = Value at lambda = forwards (also for zwanzig singlecopy)
c
c     (C) Christopher Woods     28 January 2002
c
c##############################################################

      integer ipar0,ipar1
      character*2 amber
      integer protonNum
      integer i
      double precision charge,sigma,epslon
      double precision perturbForceField
      double precision chg0,chg1,sig0,sig1,eps0,eps1
      double precision surfradius,surftension
      double precision rad0,rad1,tens0,tens1
      double precision getgbiborn,getgbscale
      double precision born0,born1,scale0,scale1

c     return 0 if something went wrong!
      getCLJ = 0
      if (ipar0.eq.ipar1) then
        getCLJ = ipar0
        return
      endif

      if (ipar0.lt.0 .or. ipar1.lt.0) then
        write(printstring,*) "Problem in setvariableparam | ",ipar0," ",ipar1
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is probably a programming bug!")
        return
      endif

c     we need to look up to see if this perturbation pair has been uused before
      do i=1,NVariableCLJ
        
        if (VariableCLJ(i,1) .eq. ipar0 .and.
     .      VariableCLJ(i,2) .eq. ipar1) then

c         we have already processed this combination
          getCLJ = -i
          return
        endif
      enddo
      
c     we need to create a new variable parameter
      NVariableCLJ = NVariableCLJ + 1
      if (NVariableCLJ.gt.MAXVARIABLEPARAM) then
        call printLine(WARNING,"Too many variable clj parameters! Increase MAXVARIABLEPARAM and recompile!")
        NVariableCLJ = NVariableCLJ - 1
        return
      endif
      
      VariableCLJ(NVariableCLJ,1) = ipar0
      VariableCLJ(NVariableCLJ,2) = ipar1
      
c     generate the parameters for these values...
      chg0 = charge(ipar0)
      chg1 = charge(ipar1)
      sig0 = sigma(ipar0)
      sig1 = sigma(ipar1)
      eps0 = epslon(ipar0)
      eps1 = epslon(ipar1)
      
      VariableChg(NVariableCLJ,1) = chg0
      VariableChg(NVariableCLJ,3) = chg1
      VariableChg(NVariableCLJ,2) = perturbForceField(chg0,chg1,lambda)

      VariableSig(NVariableCLJ,1) = sig0
      VariableSig(NVariableCLJ,3) = sig1
      VariableSig(NVariableCLJ,2) = perturbForceField(sig0,sig1,lambda)

      VariableEps(NVariableCLJ,1) = eps0
      VariableEps(NVariableCLJ,3) = eps1
      VariableEps(NVariableCLJ,2) = perturbForceField(eps0,eps1,lambda)

      VariableAtm(NVariableCLJ,1) = amber(ipar0)
      VariableAtm(NVariableCLJ,2) = amber(ipar1)

      VariableProtNum(NVariableCLJ,1) = protonNum(ipar0)
      VariableProtNum(NVariableCLJ,2) = protonNum(ipar1)

c     now generate the forwards and backwards parameters if we are doing
c     an FEP calculation (4 is backwards, 5 is forwards)
      VariableChg(NVariableCLJ,4) = perturbForceField(chg0,chg1,lambdaB)
      VariableChg(NVariableCLJ,5) = perturbForceField(chg0,chg1,lambdaF)
      VariableEps(NVariableCLJ,4) = perturbForceField(eps0,eps1,lambdaB)
      VariableEps(NVariableCLJ,5) = perturbForceField(eps0,eps1,lambdaF)
      VariableSig(NVariableCLJ,4) = perturbForceField(sig0,sig1,lambdaB)
      VariableSig(NVariableCLJ,5) = perturbForceField(sig0,sig1,lambdaF)
      
      
      if (doSA) then 
        rad0 = surfradius(ipar0)
        rad1 = surfradius(ipar1)
        tens0 = surftension(ipar0)
        tens1 = surftension(ipar1)
        
        VariableSurfRadius(NVariableCLJ,1) = rad0
        VariableSurfRadius(NVariableCLJ,3) = rad1
        VariableSurfRadius(NVariableCLJ,2) = perturbForceField(rad0,rad1,lambda)
        VariableSurfRadius(NVariableCLJ,4) = perturbForceField(rad0,rad1,lambdaB)
        VariableSurfRadius(NVariableCLJ,5) = perturbForceField(rad0,rad1,lambdaF)
        
        VariableSurfTension(NVariableCLJ,1) = tens0
        VariableSurfTension(NVariableCLJ,3) = tens1
        VariableSurfTension(NVariableCLJ,2) = perturbForceField(tens0,tens1,lambda)
        VariableSurfTension(NVariableCLJ,4) = perturbForceField(tens0,tens1,lambdaB)
        VariableSurfTension(NVariableCLJ,5) = perturbForceField(tens0,tens1,lambdaF)
        
      endif

      if (doGB) then 
        born0 = getgbiborn(ipar0)
        born1 = getgbiborn(ipar1)
        scale0 = getgbscale(ipar0)
        scale1 = getgbscale(ipar1)
        
        VariableIntrBorn(NVariableCLJ,1) = born0
        VariableIntrBorn(NVariableCLJ,3) = born1
        VariableIntrBorn(NVariableCLJ,2) = perturbForceField(born0,born1,lambda)
        
        VariablePDAScale(NVariableCLJ,1) = scale0
        VariablePDAScale(NVariableCLJ,3) = scale1
        VariablePDAScale(NVariableCLJ,2) = perturbForceField(scale0,scale1,lambda)
        VariableIntrBorn(NVariableCLJ,4) = perturbForceField(born0,born1,lambdaB)
        VariableIntrBorn(NVariableCLJ,5) = perturbForceField(born0,born1,lambdaF)
        VariablePDAScale(NVariableCLJ,4) = perturbForceField(scale0,scale1,lambdaB)
        VariablePDAScale(NVariableCLJ,5) = perturbForceField(scale0,scale1,lambdaF)
        
      endif

c     return the negative of the index...
      getCLJ = -NVariableCLJ
      
      return
      
      end
