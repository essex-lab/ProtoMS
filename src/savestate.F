

      subroutine saveState()
      
      include 'dimensions.inc'
      include 'enums.inc'
      include 'pressure.inc'
      include 'flags.inc'
      include 'printer.inc'
      include 'moves.inc'
      include 'stacks.inc'
      include 'simparams.inc'
      include 'energies.inc'
      include 'metropolis.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'tempstate.inc'
c################################################################
c
c  This routine saves the coordinates of the system (the state)
c  in temporary arrays
c
c  Samuel Genheden - 28 September 2013
c
c################################################################
      integer i,j,k
      tGCSCavity = GCSCavity
      if (boundary.eq.PERIODICBOUNDARIES) then
        do i=1,3
          tBoxDimension(i) = BoxDimension(i)
          tBoxOrigin(i) = BoxOrigin(i)
        enddo
      endif
      
c     now write out the coordinates of each protein
      do i=1,NProteins     
        do j=1,NZAtoms(i)
          do k=1,3
            tZCoords(i,j,k) = ZCoords(i,j,k)
          enddo
        enddo
      enddo

c     now write the coordinates of all of the solutes
      do i=1,NSolutes
        do j=1,3
          do k=1,3
            tSolDummy(i,j,k) = SolDummy(i,j,k)
          enddo
          tSolGeomCentre(i,j) = SolGeomCentre(i,j)
        enddo
        do j=1,NSolAtoms(i)
          do k=1,3
            tZSolCoords(i,j,k) = ZSolCoords(i,j,k)
          enddo
        enddo
      enddo

c     now write the coordinates of all of the solvents
      if (NSolvents.gt.0) then
        do i=1,NSolvents
          do j=1,3
            tSvnCentre(i,j) = SvnCentre(i,j)
          enddo
          do j=1,NSvnAtoms(i)
            do k=1,3
              tZSvnCoords(i,j,k) = ZSvnCoords(i,j,k)
            enddo
          enddo
        enddo
      endif
      
c       now write the coordinates of all of the GC-solutes
      if (NGCSolutes.gt.0) then
        do i=1,NGCSolutes
         tthetaparam(i) = thetaparam(i)
          do j=1,3
            tGCSCentre(i,j) = GCSCentre(i,j)
          enddo       
          do j=1,NGCSAtoms(i)
            do k=1,3
              tZGCSCoords(i,j,k) = ZGCSCoords(i,j,k)
            enddo
          enddo
        enddo
      endif
      
      end subroutine