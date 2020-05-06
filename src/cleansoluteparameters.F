      subroutine cleanSoluteParameters(isol)
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'solutes.inc'
      include 'constants.inc'
c#############################################################
c
c     This routine cleans up the solute parameters, and removes
c     all zero bonds/angles/dihedrals, and sets the
c     parameters for any zero chg/lj pairs to 0
c
c     (C) Christopher Woods   1 March 2002
c
c#############################################################

      integer inbnds,inangs,inubys,isol,i,ipar
            
      inbnds = NSolBonds(isol)
      if (inbnds.gt.0) then
        do i=1,inbnds
          ipar = SolBond(isol,i,3)
          if (ipar.gt.0) then
            if (BndParam(ipar,1).eq.ZERO) SolBond(isol,i,3) = 0
          endif
        enddo
      endif
      
      inangs = NSolAngles(isol)
      if (inangs.gt.0) then
        do i=1,inangs
          ipar = SolAngle(isol,i,4)
          if (ipar.gt.0) then
            if (AngParam(ipar,1).eq.ZERO) SolAngle(isol,i,4) = 0
          endif
        enddo
      endif
      
c     don't check the dihedrals as they are more trouble then they're worth!
      
      inubys = NSolUreyBradleys(isol)
      if (inubys.gt.0) then
        do i=1,inubys
          ipar = SolUreyBradley(isol,i,3)
          if (ipar.gt.0) then
            if (UByParam(ipar,1).eq.ZERO) SolUreyBradley(isol,i,3) = 0
          endif
        enddo
      endif 
      
      return
      end
