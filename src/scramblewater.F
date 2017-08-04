      subroutine scrambleWater(igcs)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'solutes.inc'
      include 'enums.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'flags.inc'

c###################################################
c
c     This routine generates a novel water configuration
c     within the gcmc system region
c
c     (C) Hannah Bruce-Macdonald, Chris Cave-Ayland, May 2016
c     
c###################################################

      integer igcs,i,j,nats
      double precision ranRange
      double precision ranval,ranvec(3),newCentre(3),transVec(3)
      double precision distance
c     pick a new location in the gcmc region
      newCentre(1) = ranRange ( coordx, coordx+gridx )
      newCentre(2) = ranRange ( coordy, coordy+gridy )
      newCentre(3) = ranRange ( coordz, coordz+gridz )
    
c     get translation vector to move the GCS to its new position
      do i=1,3
        transVec(i) = newCentre(i) - GCSCentre(igcs,i)
      enddo 

c     apply trans vec to all atoms of the GCS and update its centre
      do j=1,3
        do i=1,NGCSAtoms(igcs)
          ZGCSCoords(igcs,i,j) = ZGCSCoords(igcs,i,j) + transVec(j)
        enddo
        GCSCentre(igcs,j) = newCentre(j)
      enddo

      if ( sphericalgcmc ) then

22    newCentre(1) = ranRange ( centerx-gcradius, centerx+gcradius )
      newCentre(2) = ranRange ( centery-gcradius, centery+gcradius )
      newCentre(3) = ranRange ( centerz-gcradius, centerz+gcradius )
    
c     get translation vector to move the GCS to its new position
      do i=1,3
        transVec(i) = newCentre(i) - GCSCentre(igcs,i)
      enddo 

c     apply trans vec to all atoms of the GCS and update its centre
      do j=1,3
        do i=1,NGCSAtoms(igcs)
          ZGCSCoords(igcs,i,j) = ZGCSCoords(igcs,i,j) + transVec(j)
        enddo
        GCSCentre(igcs,j) = newCentre(j)
      enddo

      distance =((GCSCentre(igcs,1)-centerx)**2+
     : (GCSCentre(igcs,2)-centery)**2+(GCSCentre(igcs,3)-centerz)**2)**0.5
      if (distance.gt.gcradius) then
      goto 22
      endif 
      endif

c     pick a randomly oriented vector and rotate around it in the range -pi to pi
c     completely scrambles the orientation
      call getVectorOnSphere(ranvec)
      ranval = ranRange( -PI, PI ) 
      call rotateGCSolute(igcs,ranvec,ranval)
      
      return
      end
