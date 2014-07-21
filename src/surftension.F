      double precision function surftension(ipar)
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'constants.inc'
      include 'surface.inc'
c###############################################
c
c     Return surface tension of simulation parameter
c     ipar at the current value of lambda
c
c     Julien Michel November 2004
c
c###############################################
      integer ipar
      surftension  = ZERO
      
      if (ipar.eq.0 .or. ipar.gt.NCLJTypes) then
        return
      else if (ipar.gt.0) then
        surftension = SurfTensionDB(CLJType(ipar))
      else if (-ipar .le. NVariableCLJ) then
        surftension = VariableSurfTension(-ipar,2)
      endif
      
      return
      end
