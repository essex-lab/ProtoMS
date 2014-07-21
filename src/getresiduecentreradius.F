      subroutine getResidueCentreRadius(iprot,ires)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'proteins.inc'
c#################################################################
c
c     This routine will calculate the centre and radius
c     of residue ires in protein iprot
c
c     (C) Christopher Woods   6 March 2002
c
c#################################################################

      integer iprot,ires
      integer i,istrt,iend,j
      double precision imax(3),imin(3),dist2
      
      if (iprot.le.0 .or. iprot.gt.NProteins) return
      if (ires.le.0 .or. ires.gt.NResidues(iprot)) return
      
c     assume residue runs concurrently from the nitrogen to the oxygen
      istrt = BBZid(iprot,ires,1)
      iend = BBZid(iprot,ires,4)

      do i=1,3
        imax(i) = -999999.0
        imin(i) = 999999.0
      enddo

      resRadius(iprot,ires) = ZERO

c     calculate the centre      
      do i=istrt,iend
c       don't include dummy atoms
        if (ZPar(iprot,i).ne.0) then
          do j=1,3
            if (ZCoords(iprot,i,j).gt.imax(j)) imax(j) = ZCoords(iprot,i,j)
            if (ZCoords(iprot,i,j).lt.imin(j)) imin(j) = ZCoords(iprot,i,j)
          enddo
        endif
      enddo

c     the centre is halfway between the max and min values
      do i=1,3
        ResCentre(iprot,ires,i) = ((imax(i)-imin(i)) * 0.5D+00) + imin(i)
      enddo

c     the radius(squared) is the maximum distance(squared) between 
c     the centre and an atom
      do i=istrt,iend
c       don't include dummy atoms
        if (ZPar(iprot,i).ne.0) then
          dist2 = (ZCoords(iprot,i,1)-ResCentre(iprot,ires,1))**2 +
     .            (ZCoords(iprot,i,2)-ResCentre(iprot,ires,2))**2 +
     .            (ZCoords(iprot,i,3)-ResCentre(iprot,ires,3))**2
          if (dist2.gt.resRadius(iprot,ires)) resRadius(iprot,ires) = dist2
        endif
      enddo

      resRadius(iprot,ires) = dsqrt(resRadius(iprot,ires))

      return
      end
