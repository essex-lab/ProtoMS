      double precision function minSolventProteinDist2(iprot,isvn)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'enums.inc'
      include 'simparams.inc'
      include 'flags.inc'
      include 'proteins.inc'
c##############################################################
c
c     This routine estimates the minimum distance squared
c     between solvent isvn and protein iprot (it will always
c     overestimate, never underestimate)
c
c     (C) Christopher Woods
c
c##############################################################

      integer isvn,iprot
      double precision distx,disty,distz,dist2
      integer i
      
      minSolventProteinDist2 = 1.D+10

      if (NResidues(iprot).le.0) return
      
      do i=1,NResidues(iprot)
        distx = dabs(ResCentre(iprot,i,1)-SvnCentre(isvn,1))
        disty = dabs(ResCentre(iprot,i,2)-SvnCentre(isvn,2))
        distz = dabs(ResCentre(iprot,i,3)-SvnCentre(isvn,3))
        if (boundary.eq.PERIODICBOUNDARIES) then
          if (distx.gt.BoxDimension2(1)) then
            distx = distx - BoxDimension(1)
          endif
          if (disty.gt.BoxDimension2(2)) then
            disty = disty - BoxDimension(2)
          endif
          if (distz.gt.BoxDimension2(3)) then
            distz = distz - BoxDimension(3)
          endif
        endif
        dist2 = distx**2 + disty**2 + distz**2
        if (dist2.lt.minSolventProteinDist2) 
     .       minSolventProteinDist2 = dist2
      enddo
      
      return
      end
