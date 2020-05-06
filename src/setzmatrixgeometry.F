
      subroutine setZmatrixGeometry(ityp,imol,ires,igeom)
      include 'dimensions.inc'
      include 'printer.inc'
c##############################################################
c
c     This routine will eventually do exactly the same
c     as setEquilibriumGeometry, but it will reset
c     degrees of freedom back to their z-matrix defined
c     values
c
c     (C) Christopher Woods, January 2005
c
c##############################################################

      integer ityp,imol,ires,igeom
      
      call printLine(WARNING,"setZmatrixGeometry has not yet been implemented!")
      
      return
      end
