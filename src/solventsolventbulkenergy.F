      
      double precision function solventSolventBulkEnergy(cnrg,ljnrg,
     .                        cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'constants.inc'
      include 'flags.inc'
c##################################################################
c
c     This function returns the total solvent-solvent energy, together
c     with the coulomb and lj components in cnrg and ljnrg.
c
c     (C) Christopher Woods    29 January 2002
c
c##################################################################

      double precision cnrg,ljnrg,nrg,icnrg,iljnrg,inrg
      double precision cnrgf,ljnrgf,cnrgb,ljnrgb,icnrgf,icnrgb,iljnrgf,iljnrgb
      double precision mindist,imindist
      double precision solventSolventSingleEnergy
      integer i,j

      solventSolventBulkEnergy = ZERO
      mindist = 99990.0
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      nrg = ZERO

      if (NSolvents.lt.2) return
      
      do i=1,NSolvents-1
        if (.not.DoneSvnNrg(i)) then
          do j=i+1,NSolvents
            if (.not.DoneSvnNrg(j)) then
              inrg = solventSolventSingleEnergy(i,j,
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

      solventSolventBulkEnergy = nrg
      
      return
      
      end
