      subroutine rotateGCSolute(igcs,ranvec,ranval)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'solutes.inc'
      include 'constants.inc'
      integer igcs,i,j
      double precision ranvec(3),ranval,ivec(3)
      double precision vecLength,ilgth
      integer nats

      nats = NGCSAtoms(igcs)
      do j=1,nats
        do i=1,3
          ivec(i) = ZGCSCoords(igcs,j,i) - GCSCentre(igcs,i)
        enddo
        ilgth = vecLength(ivec)
        if (ilgth .gt. ZERO) then
          call rotateVector(ivec,ranvec,ranval)
c         restore the length of the vector - prevents precision
c         errors from propagating
          call setLength(ilgth,ivec)
          do i=1,3
            ZGCSCoords(igcs,j,i) = GCSCentre(igcs,i) + ivec(i)
          enddo
        endif
      enddo
            
      return
      end
