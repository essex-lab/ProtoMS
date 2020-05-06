      double precision function minSolventSoluteDist2(isvn,isol)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'enums.inc'
      include 'simparams.inc'
      include 'solutes.inc'
      include 'flags.inc'
c##############################################################
c
c     This routine estimates the minimum distance squared
c     between solvent isvn and solute isol (it will always
c     overestimate, never underestimate)
c
c     (C) Christopher Woods
c
c##############################################################

      integer isvn,isol
      double precision distx,disty,distz
      
      minSolventSoluteDist2 = 1.D+10
      
      distx = dabs(SolCentre(isol,1)-SvnCentre(isvn,1))
      disty = dabs(SolCentre(isol,2)-SvnCentre(isvn,2)) 
      distz = dabs(SolCentre(isol,3)-SvnCentre(isvn,3))
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
      
      minSolventSoluteDist2 = distx**2 + disty**2 + distz**2

      return
      end
