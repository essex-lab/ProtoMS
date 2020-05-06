      
      double precision function GCSoluteBulkEnergy(cnrg,ljnrg,
     .                        cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'constants.inc'
      include 'flags.inc'
c##################################################################
c
c     This function returns the total GCS energy, together
c     with the coulomb and lj components in cnrg and ljnrg.
c
c     (C) MSB Sepy 2012
c
c##################################################################

      double precision cnrg,ljnrg,nrg,icnrg,iljnrg,inrg
      double precision cnrgf,ljnrgf,cnrgb,ljnrgb,icnrgf,icnrgb,iljnrgf,iljnrgb
      double precision mindist,imindist
      double precision GCSoluteGCSoluteSingleEnergy
      integer i,j

      GCSoluteBulkEnergy = ZERO
      mindist = 99990.0
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      nrg = ZERO

      if (NGCSolutes.lt.2) return
      
      do i=1,NGCSolutes-1
        if (.not.DoneGCSNrg(i)) then
          do j=i+1,NGCSolutes
            if (.not.DoneGCSNrg(j)) then
              inrg = GCSoluteGCSoluteSingleEnergy(i,j,
     .            icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
              nrg = nrg + inrg
              cnrg = cnrg + icnrg
              ljnrg = ljnrg + iljnrg
              cnrgf = cnrgf + icnrgf
              cnrgb = cnrgb + icnrgb
              ljnrgf = ljnrgf + iljnrgf
              ljnrgb = ljnrgb + iljnrgb
              if (imindist.lt.mindist) mindist = imindist
            endif
          enddo
        endif
      enddo

      GCSoluteBulkEnergy = nrg
      
      return
      
      end
