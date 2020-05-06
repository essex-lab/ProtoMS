      double precision function surfradius(ipar)
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'constants.inc'
      include 'surface.inc'
      include 'printer.inc'
c###############################################
c
c     Return surface radius of simulation parameter
c     ipar at the current value of lambda
c
c
c     Julien Michel November 2004
c
c###############################################
      integer ipar
      surfradius = ZERO
      
      if (ipar.eq.0 .or. ipar.gt.NCLJTypes) then
        return
      else if (ipar.gt.0) then
        surfradius = SurfRadDB(CLJType(ipar))
      else if (-ipar .le. NVariableCLJ) then
        surfradius = VariableSurfRadius(-ipar,2)
      endif
      
      return
      end
