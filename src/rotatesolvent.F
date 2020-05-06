      subroutine rotateSolvent(isvn,ranvec,ranval)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'constants.inc'
      integer isvn,i,j
      double precision ranvec(3),ranval,ivec(3)
      double precision vecLength,ilgth
      integer nats

      nats = NSvnAtoms(isvn)
      do j=1,nats
        do i=1,3
          ivec(i) = ZSvnCoords(isvn,j,i) - SvnCentre(isvn,i)
        enddo
        ilgth = vecLength(ivec)
        if (ilgth .gt. ZERO) then
          call rotateVector(ivec,ranvec,ranval)
c         restore the length of the vector - prevents precision
c         errors from propagating
          call setLength(ilgth,ivec)
          do i=1,3
            ZSvnCoords(isvn,j,i) = SvnCentre(isvn,i) + ivec(i)
          enddo
        endif
      enddo
            
      return
      end
