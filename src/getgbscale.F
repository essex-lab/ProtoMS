      double precision function getgbscale(ipar)
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'constants.inc'
      include 'gb.inc'
      include 'printer.inc'
c###############################################
c
c     Return PDA scale factor of simulation parameter
c     ipar at the current value of lambda
c
c
c     Julien Michel November 2004
c
c###############################################
      integer ipar
      getgbscale = ZERO
      
      if (ipar.eq.0 .or. ipar.gt.NCLJTypes) then
        return
      else if (ipar.gt.0) then
        getgbscale = PDAScaleDB(CLJType(ipar))
      else if (-ipar .le. NVariableCLJ) then
        getgbscale = VariablePDAScale(-ipar,2)
      endif
      
      return
      end
