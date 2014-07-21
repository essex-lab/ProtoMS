      subroutine getSoluteCentreRadius(isol)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'solutes.inc'
c#############################################################
c
c     This routine calculates the centre and radius of
c     solute isol
c
c     (C) Christopher Woods  6 March 2002
c
c#############################################################

      integer isol,i,nats,j
      double precision imax(3),imin(3),dist2
      
      if (isol.le.0 .or. isol.gt.NSolutes) return
      nats = NSolAtoms(isol)
      if (nats.le.0) return
      
      do i=1,3
        imin(i) = 9999900.0
        imax(i) = -999999.0
      enddo

      solRadius(isol) = ZERO

c     get the maximum and minimum x/y/z values
      do i=1,nats
c       don't include dummy atoms
        if (ZSolPar(isol,i).ne.0) then
          do j=1,3
            if (ZSolCoords(isol,i,j).gt.imax(j)) imax(j) = ZSolCoords(isol,i,j)
            if (ZSolCoords(isol,i,j).lt.imin(j)) imin(j) = ZSolCoords(isol,i,j)
          enddo
        endif
      enddo

c     the centre is halfway between the max and min values
      do j=1,3
        SolCentre(isol,j) = ((imax(j)-imin(j)) * 0.5D+00 + imin(j))
      enddo

c     the radius is the largest distance of an atom from the centre
      do i=1,nats
c       don't include dummy atoms
        if (ZSolPar(isol,i).ne.0) then
          dist2 = (ZSolCoords(isol,i,1) - SolCentre(isol,1))**2 + 
     .            (ZSolCoords(isol,i,2) - SolCentre(isol,2))**2 + 
     .            (ZSolCoords(isol,i,3) - SolCentre(isol,3))**2
          if (dist2.gt.SolRadius(isol)) SolRadius(isol) = dist2
        endif
      enddo

      SolRadius(isol) = dsqrt(SolRadius(isol))

      return
      end
