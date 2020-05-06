      double precision function getgbiborn(ipar)
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'constants.inc'
      include 'gb.inc'
      include 'printer.inc'
c###############################################
c
c     Return intrinsic born radius of simulation parameter
c     ipar at the current value of lambda
c
c
c     Julien Michel November 2004
c
c###############################################
      integer ipar
      getgbiborn = ZERO
      
      if (ipar.eq.0 .or. ipar.gt.NCLJTypes) then
        return
      else if (ipar.gt.0) then
        getgbiborn = IntrBornDB(CLJType(ipar))
      else if (-ipar .le. NVariableCLJ) then
        getgbiborn = VariableIntrBorn(-ipar,2)
      endif
      
      return
      end
