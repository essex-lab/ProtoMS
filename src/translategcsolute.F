      subroutine translateGCSolute(igcs,idx)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'solutes.inc'
      include 'simparams.inc'
      integer igcs,i,j
      double precision idx(3)
      integer nats

1      nats = NGCSAtoms(igcs)
      do j=1,3
        do i=1,nats
          ZGCSCoords(igcs,i,j) = ZGCSCoords(igcs,i,j) + idx(j) 
        enddo 
          GCSCentre(igcs,j) = GCSCentre(igcs,j) + idx(j) 
        
      enddo

      return
      end
