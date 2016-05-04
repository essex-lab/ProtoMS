      subroutine scrambleWater(igcs,ieold,ienew)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'solutes.inc'
      include 'enums.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'flags.inc'
c###################################################
c
c     This routine moves solvent molecule igcs
c     via translation and rotation
c
c     Translation and rotation is done in the 
c     within the GC cavity
c
c     (C) MSB May 09
c     
c     
c###################################################

      integer igcs,i,j,nats
      double precision ranRange,ienew,ieold
      double precision ranval,ranvec(3),newCentre(3),transVec(3)
      double precision pivalue     
      pivalue = 4 * atan( 1.0 )
      ienew = 0

      newCentre(1) =ranRange ( coordx, coordx+gridx )
      newCentre(2) =ranRange ( coordy, coordy+gridy )
      newCentre(3) =ranRange ( coordz, coordz+gridz )
    
 
      do i=1,3
        transVec(i) = newCentre(i) - GCSCentre(igcs,i)
      enddo 

      nats = NGCSAtoms(igcs)
      do j=1,3
        do i=1,nats
          ZGCSCoords(igcs,i,j) = ZGCSCoords(igcs,i,j) + transVec(j)
        enddo
        GCSCentre(igcs,j) = newCentre(j)
      enddo


c     translate the GCsolute molecule randomly
      call getVectorOnSphere(ranvec)
c      ranval = ranRange(-GCSADel(igcs),GCSADel(igcs))
      
      ranval = ranRange(- pivalue , pivalue ) 
     
 
c      do i=1,3
c        ranvec(i) = ranval * ranvec(i)
c      enddo
c      call translateGCSolute(igcs,ranVec)
c      PRINT*, "Post moved" ,GCSCentre(igcs,1) ,GCSCentre(igcs,2) ,GCSCentre(igcs,3)
c     now rotate the GCsolute
c      call getVectorOnSphere(ranvec)
c     ranval = ranRange(-GCSADel(igcs),GCSADel(igcs))

      call rotateGCSolute(igcs,ranvec,ranval)

c       PRINT*, "Post moved" ,GCSCentre(igcs,1) ,GCSCentre(igcs,2) ,GCSCentre(igcs,3)

c        if (GCSCentre(igcs,1).gt.(coordx+gridx+TINY)) then
c        ienew = (ienew + (HUGE*HUGE))
c        endif
c        if (GCSCentre(igcs,2).gt.(coordy+gridy+TINY)) then
c        ienew = (ienew + (HUGE*HUGE))
c         endif
c        if (GCSCentre(igcs,3).gt.(coordz+gridz+TINY)) then
c        ienew = (ienew + (HUGE*HUGE))
c        endif  
c        if (GCSCentre(igcs,1).lt.(coordx-TINY)) then
c        ienew = (ienew + (HUGE*HUGE))
c         endif
c       if (GCSCentre(igcs,2).lt.(coordy-TINY)) then
c        ienew = (ienew + (HUGE*HUGE))
c        endif
c        if (GCSCentre(igcs,3).lt.(coordz-TINY)) then
c        ienew = (ienew + (HUGE*HUGE))
c        endif
      return
      end
