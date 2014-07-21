      integer function cljID(ipar)
      include 'dimensions.inc'
      include 'simparams.inc'
      include 'parameters.inc'
      include 'constants.inc'
c#############################################
c
c     This function returns the index in the
c     CLJ database of the simulation clj
c     parameter zpar, or 0 if it doesn't
c     exist or is a dummy - if this is 
c     a variable parameter then it returns
c     the parameter at lam=0 for lam<=0.5,
c     and the parameter at lam=1 for lam>0.5
c
c     (C) Christopher Woods
c
c#############################################

      integer ipar
            
      if (ipar.eq.0 .or. ipar.gt.NCLJTypes) then
        cljID = 0
      else if (ipar.lt.0) then
        if (-ipar .gt. NVariableCLJ) then
          cljID = 0
        else if (lambda.le.HALF) then
          if (VariableCLJ(-ipar,1).eq.0) then
            cljID = 0
          else
            cljID = CLJType(VariableCLJ(-ipar,1))
          endif
        else
          if (VariableCLJ(-ipar,2).eq.0) then
            cljID = 0
          else
            cljID = CLJType(VariableCLJ(-ipar,2))
          endif
        endif
      else      
        cljID = CLJType(ipar)
      endif
        
      return
      end
