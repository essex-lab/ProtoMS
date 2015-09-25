      double precision function solventSolventEnergy(savecut,isvn,
     .             cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'constants.inc'
      include 'flags.inc'
c##############################################################
c
c     This function returns the interaction energy between 
c     solvent molecule isvn, and all other solvent molecules...
c
c     This function will cut the interaction if the distance
c     between the solvent centres is greater than cut
c
c     (C) Christopher Woods     29 January 2002
c
c##############################################################

      integer isvn,savecut
      double precision mindist,imindist
      double precision cnrg,ljnrg,icnrg,iljnrg,nrg,inrg
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb,icnrgf,icnrgb,iljnrgf,iljnrgb
      double precision solventSolventSingleEnergy
      integer i,nats
      
      solventSolventEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      nrg = ZERO
      mindist = 999990.0
      
      if (NSolvents.lt.2 .or. isvn.gt.NSolvents) return
      nats = NSvnAtoms(isvn)
      if (nats.le.0) return
      
      if (DoneSvnNrg(isvn)) return
      
!$OMP PARALLEL DO DEFAULT(NONE)

!$OMP&PRIVATE(inrg,icnrg,iljnrg,icnrgf,icnrgb,iljnrgf,iljnrgb,imindist)

!$OMP&SHARED(isvn)

!$OMP&REDUCTION(+:nrg,cnrg,ljnrg,cnrgf,cnrgb,ljnrgf,ljnrgb)
!$OMP&REDUCTION(min:mindist)
      do i=1,NSolvents
        if (i.ne.isvn) then
          inrg =  solventSolventSingleEnergy(isvn,i,
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

      solventSolventEnergy = nrg
      return
      end
