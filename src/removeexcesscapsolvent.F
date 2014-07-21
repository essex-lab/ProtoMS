
      subroutine removeExcessCapSolvent(ox,oy,oz,irad)
      include 'dimensions.inc'
      include 'solvents.inc'
c##########################################################
c
c     This routine removes any solvent molecules
c     that have their centers outside the sphere
c     centered at (ox,oy,oz), radius irad.
c
c     (C) Christopher Woods, Febrauary 2005
c
c##########################################################

      double precision ox,oy,oz,irad
      double precision irad2,dist2
      integer i,nsvn
      

      nsvn = NSolvents
      irad2 = irad**2
      do i=nsvn,1,-1
        dist2 = (SvnCentre(i,1) - ox)**2 + (SvnCentre(i,2) -oy)**2 +
     .    (SvnCentre(i,3)-oz)**2
        if (dist2 .gt. irad2) then 
          goto 100
        else
          goto 110
        endif
 100    continue
c       delete this solvent molecule        
        if (i.eq.NSolvents) then
          NSolvents = NSolvents - 1
        else
          call copySolvent(NSolvents,i)
          NSolvents = NSolvents - 1
        endif
 110    continue
c       don't delete the solvent molecule        
      enddo
      
      return
      end
      
