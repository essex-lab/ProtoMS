      subroutine moveGCSolute(igcs,ieold,ienew)
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

      integer igcs,i
      double precision ranRange,ienew,ieold
      double precision ranval,ranvec(3)
      double precision distance

      ienew = 0
c     translate the GCsolute molecule randomly
      call getVectorOnSphere(ranvec)
      ranval = ranRange(-GCSADel(igcs),GCSADel(igcs))
      do i=1,3
        ranvec(i) = ranval * ranvec(i)
      enddo
      call translateGCSolute(igcs,ranVec)
c      PRINT*, "Post moved" ,GCSCentre(igcs,1) ,GCSCentre(igcs,2) ,GCSCentre(igcs,3)
c     now rotate the GCsolute
      call getVectorOnSphere(ranvec)
      ranval = ranRange(-GCSRDel(igcs),GCSRDel(igcs))
      call rotateGCSolute(igcs,ranvec,ranval)
c       PRINT*, "Post moved" ,GCSCentre(igcs,1) ,GCSCentre(igcs,2) ,GCSCentre(igcs,3)
c calculate the distance

        if ( sphericalgcmc ) then
        print*, 'movegcsolute'
        distance = ((GCSCentre(igcs,1)-centerx)**2.+(GCSCentre(igcs,2)
     :  -centery)**2+(GCSCentre(igcs,3)-centerz)**2)**0.5
        if (distance.gt.gcradius) then
        ienew = (ienew + (HUGE*HUGE))
        endif
        
        else
        if (GCSCentre(igcs,1).gt.(coordx+gridx+TINY)) then
        ienew = (ienew + (HUGE*HUGE))
        endif
        if (GCSCentre(igcs,2).gt.(coordy+gridy+TINY)) then
        ienew = (ienew + (HUGE*HUGE))
         endif
        if (GCSCentre(igcs,3).gt.(coordz+gridz+TINY)) then
        ienew = (ienew + (HUGE*HUGE))
        endif  
        if (GCSCentre(igcs,1).lt.(coordx-TINY)) then
        ienew = (ienew + (HUGE*HUGE))
         endif
       if (GCSCentre(igcs,2).lt.(coordy-TINY)) then
        ienew = (ienew + (HUGE*HUGE))
        endif
        if (GCSCentre(igcs,3).lt.(coordz-TINY)) then
        ienew = (ienew + (HUGE*HUGE))
        endif
        endif
      return
      end
