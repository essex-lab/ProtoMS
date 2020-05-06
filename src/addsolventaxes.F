      subroutine addSolventAxes(isvn)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'simparams.inc'
      include 'solvents.inc'
c########################################################
c     
c     This routine gets the major, minor and perpendicular
c     axes of rotation, and the mean centre of solvent
c     isvn - this is useful for performing solvent
c     translations and rotations
c
c     (C) Christopher Woods  4 March 2002
c
c########################################################

      integer isvn,i,j
      double precision itot(3),imaj(3),imin(3),iper(3)
      double precision max,min,coords
      
      if (isvn.le.0 .or. isvn.gt.NSolvents) return

      if ( .not. SIRE_COMPATIBILITY_MODE ) then
        i = isvn
        call getMajorMinorAxis(SOLVENT,i,0,itot,imaj,imin,iper)
      
        do 100 i=1,3
          SvnCentre(isvn,i) = itot(i)
c          SvnMajAxis(isvn,i) = imaj(i)
c          SvnMinAxis(isvn,i) = imin(i)
c          SvnPerAxis(isvn,i) = iper(i)
100     continue

      else

        do 200 i=1,3
          max = -99999999.0
          min = 99999999.0

          do 201 j=1,NSvnAtoms(isvn)
            coords = ZSvnCoords(isvn,j,i)

            if (coords .lt. min) min = coords
            if (coords .gt. max) max = coords
201       continue

          SvnCentre(isvn,i) = min + 0.5*(max-min)
200     continue

      endif

      return
      end

