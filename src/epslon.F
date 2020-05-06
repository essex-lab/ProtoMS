      double precision function epslon(ipar)
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'constants.inc'
c###############################################
c
c     Return epsilon of simulation parameter
c     ipar at the current value of lambda
c
c     (C) Christopher Woods, October 2004
c
c###############################################

      integer ipar
      
      epslon = ZERO
      
      if (ipar.eq.0 .or. ipar.gt.NCLJTypes) then
        return
      else if (ipar.gt.0) then
        epslon = EpsDB(CLJType(ipar))
      else if (-ipar .le. NVariableCLJ) then
        epslon = VariableEps(-ipar,2)
      endif
      
      return
      end
