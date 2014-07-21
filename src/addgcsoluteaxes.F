      subroutine addGCsoluteAxes(igcs)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'solvents.inc'
c########################################################
c     
c     This routine gets the major, minor and perpendicular
c     axes of rotation, and the mean centre of solute
c     igcs - this is useful for performing solute
c     translations and rotations
c
c     MSB 1st April 09
c
c########################################################

      integer igcs,i
      double precision itot(3),imaj(3),imin(3),iper(3)
      
      if (igcs.le.0 .or. igcs.gt.NGCSolutes) return
      i = igcs
      call getMajorMinorAxis(GCSOLUTE,i,0,itot,imaj,imin,iper)
      
      do 100 i=1,3
        GCSCentre(igcs,i) = itot(i)
c        SvnMajAxis(isvn,i) = imaj(i)
c        SvnMinAxis(isvn,i) = imin(i)
c        SvnPerAxis(isvn,i) = iper(i)
100   continue
c       PRINT*, "Centre-GCS" ,GCSCentre(igcs,1),GCSCentre(igcs,2),GCSCentre(igcs,3)
      return
      end
