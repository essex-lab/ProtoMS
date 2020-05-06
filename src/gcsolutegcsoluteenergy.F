      double precision function GCsoluteGCsoluteEnergy(savecut,igcs,
     .             cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'stacks.inc'
      include 'surface.inc'
      include 'constants.inc'
      include 'printer.inc'
      include 'flags.inc'
c##############################################################
c
c     This function returns the interaction energy between 
c     gc-solute molecule igcs, and all other gcs's
c
c     This function will cut the interaction if the distance
c     between the solute centres is greater than cut
c
c     (C) MSB 22/04/09
c
c##############################################################

      integer igcs,savecut
      double precision mindist,imindist
      double precision cnrg,ljnrg,icnrg,iljnrg,nrg,inrg
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb,icnrgf,icnrgb,iljnrgf,iljnrgb
      double precision GCsoluteGCsoluteSingleEnergy
      integer i,nats
      
      GCsoluteGCsoluteEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      nrg = ZERO
      mindist = 999990.0
      
      if (NGCSolutes.lt.2 .or. igcs.gt.NGCSolutes) return
                       
      nats = NGCSAtoms(igcs)
      
      if (nats.le.0) return

      if (DoneGCSNrg(igcs)) return
      
      do i=1,NGCSolutes
        if (i.ne.igcs) then
          inrg =  GCsoluteGCsoluteSingleEnergy(igcs,i,
     .                 icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
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

      GCsoluteGCsoluteEnergy = nrg
c      PRINT*, "GCS-GCS",GCsoluteGCsoluteEnergy 
      return
      end
