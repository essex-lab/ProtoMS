      double precision function proteinSoluteEnergy(savecut,iprot,isol,cnrg,ljnrg,
     .  cnrgb,ljnrgb,cnrgf,ljnrgf,mindist,scl,sclf,sclb)
      include 'dimensions.inc'
      include 'proteins.inc'
      include 'simparams.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'enums.inc'
c#############################################################
c
c     This function returns the total interaction energy between
c     protein iprot and solute isol
c
c     (C) Christopher Woods      29 January 2002
c
c#############################################################

      integer iprot,isol,savecut
      double precision cnrg,ljnrg,nrg,icnrg,iljnrg
      double precision cnrgf,cnrgb,ljnrgf,ljnrgb,icnrgf,icnrgb,iljnrgf,iljnrgb
      double precision mindist,imindist,sclfac
      integer i,nres
      double precision residueSoluteEnergy
      double precision scl,sclf,sclb
      proteinSoluteEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      nrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      mindist = 99000.0
      
      if (NProteins.le.0 .or. iprot.gt.NProteins .or. NSolutes.le.0
     .         .or. isol.gt.NSolutes) return

c      PRINT*, "I AM",isol 
      if (IgnoreSolute(isol)) return
          
      nres = NResidues(iprot)
      if (nres.le.0) return
      
      do i=1,nres
        nrg = nrg + residueSoluteEnergy(iprot,i,isol,
     .               icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist,
     .               scl,sclf,sclb)
        cnrg = cnrg + icnrg
        ljnrg = ljnrg + iljnrg
        cnrgf = cnrgf + icnrgf
        cnrgb = cnrgb + icnrgb
        ljnrgf = ljnrgf + iljnrgf
        ljnrgb = ljnrgb + iljnrgb
        if (imindist.lt.mindist) mindist = imindist
      enddo

      if (cuttype.eq.MOLECULEMOLECULE) then
c       performing protein-molecule solute-molecule cutting
        if (mindist.lt.cut2) then
c         within the cutoff distance
          if (mindist.gt.feather2) then
c           however, in the feather region
            sclfac = featherfac * (cut2-mindist)
            nrg = nrg * sclfac
            cnrg = cnrg * sclfac
            ljnrg = ljnrg * sclfac
            cnrgf = cnrgf * sclfac
            ljnrgf = ljnrgf * sclfac
            cnrgb = cnrgb * sclfac
            ljnrgb = ljnrgb * sclfac
          endif
        else
          nrg = ZERO
          cnrg = ZERO
          ljnrg = ZERO
          cnrgf = ZERO
          cnrgb = ZERO
          ljnrgf = ZERO
          ljnrgb = ZERO
        endif
      endif
      
      proteinSoluteEnergy = nrg
      
      return
      end
