      double precision function charge(ipar)
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'constants.inc'
c###############################################
c
c     Return charge of simulation parameter
c     ipar at the current value of lambda
c
c     (C) Christopher Woods, October 2004
c
c###############################################

      integer ipar
      
      charge = ZERO
      
      if (ipar.eq.0 .or. ipar.gt.NCLJTypes) then
        return
      else if (ipar.gt.0) then
        charge = ChgDB(CLJType(ipar))
      else if (-ipar .le. NVariableCLJ) then
        charge = VariableChg(-ipar,2)
      endif
      
      return
      end
