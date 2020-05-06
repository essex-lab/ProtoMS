      subroutine translateSolvent(isvn,idx)
      include 'dimensions.inc'
      include 'solvents.inc'
      integer isvn,i,j
      double precision idx(3)
      integer nats

      nats = NSvnAtoms(isvn)
      do j=1,3
        do i=1,nats
          ZSvnCoords(isvn,i,j) = ZSvnCoords(isvn,i,j) + idx(j)
        enddo 
          SvnCentre(isvn,j) = SvnCentre(isvn,j) + idx(j)
      enddo

      return
      end
