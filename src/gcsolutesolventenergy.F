      
      double precision function GCsoluteSolventEnergy(igcs,
     .            cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist,
     .            scl,sclf,sclb)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'flags.inc'
c################################################################
c
c     This function calculates the total interaction energy
c     between GCsolute igcs, and all the solvent
c
c     (C) MSb 22nd April 2009
c
c################################################################
      include 'printer.inc'
      integer igcs
      integer i
      double precision cnrg,ljnrg,nrg,icnrg,iljnrg,mindist,imindist,inrg
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb,icnrgf,icnrgb,iljnrgf,iljnrgb
      double precision GCsoluteSingleSolventEnergy
      double precision scl,sclf,sclb

      GCsoluteSolventEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      nrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      mindist = 999999.0

      if (igcs.gt.NGCSolutes  .or. NGCSolutes.le.0 .or. NSolvents.le.0) return

      do i=1,NSolvents
        inrg = GCsoluteSingleSolventEnergy(igcs,i,
     .         icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist,
     .         scl,sclf,sclb)
        nrg = nrg + inrg
        cnrg = cnrg + icnrg
        ljnrg = ljnrg + iljnrg
        cnrgf = cnrgf + icnrgf
        cnrgb = cnrgb + icnrgb
        ljnrgf = ljnrgf + iljnrgf
        ljnrgb = ljnrgb + iljnrgb 
        if (imindist.lt.mindist) mindist = imindist
      enddo

      GCsoluteSolventEnergy = nrg

      return
      end
