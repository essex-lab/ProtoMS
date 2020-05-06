      subroutine calculateSolventBox
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'flags.inc'
      include 'simparams.inc'
c##############################################
c
c     This routine generates the dimensions
c     of the solvent box (using a very basic
c     and probably wrong algorithm!)
c
c     (C) Christopher Woods, October 2004
c
c##############################################
      
      integer i,k
      
c     set the box dimension initially from the first solvent atom
      do i=1,3
        boxTop(i) = ZSvnCoords(1,1,i)
        boxOrigin(i) = ZSvnCoords(1,1,i)
      enddo

c     now find the extent of the solvent            
      do i=1,NSolvents
        do k=1,3
          if (SvnCentre(i,k).gt.boxTop(k))
     .              boxTop(k) = SvnCentre(i,k)
          if (SvnCentre(i,k).lt.boxOrigin(k))
     .              boxOrigin(k) = SvnCentre(i,k)
        enddo
      enddo
     
c     now use the origin and top of the box to obtain the sizes, and
c     (size/2)**2
      do i=1,3
        BoxDimension(i) = BoxTop(i)-BoxOrigin(i)
c       get half the box dimension (for periodic boundary test)
        BoxDimension2(i) = 0.5D+00*BoxDimension(i)
      enddo

      SolventBox = .true.

      return
      end
