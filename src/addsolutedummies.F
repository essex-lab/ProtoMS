      subroutine addSoluteDummies(isol)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'solutes.inc'
c########################################################
c
c     This adds dummy atoms to the solute.
c
c     The first is placed at the mean geometric centre
c     of the molecule, the second placed 1.0D+00 along
c     the major axis, from 1, and the third placed
c     1.0D+00 from 1 along the minor axis. The major, minor
c     and perpendicular axes are saved for use in
c     rotating the solute
c
c     (C) Christopher Woods    19 February 2002
c
c     Julien Michel October 2004
c     The major,minor,perpendicular axis are not saved anymore
c#######################################################

      integer isol
      double precision itot(3),imaj(3),imin(3),iper(3)
      integer i

      if (isol.le.0 .or. isol.gt.NSolutes) return
      
      call getMajorMinorAxis(SOLUTE,isol,0,itot,imaj,imin,iper)

      do 400 i=1,3
        SolDummy(isol,1,i) = itot(i)
        SolDummy(isol,2,i) = SolDummy(isol,1,i) + imaj(i)
        SolDummy(isol,3,i) = SolDummy(isol,1,i) + imin(i)
        SolGeomCentre(isol,i) = itot(i)
400   continue        

      return
      end
