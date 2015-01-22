      subroutine loadState()
      
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
c  This routine reads the coordinates of the system (the state)
c  from temporary arrays
c
c  Samuel Genheden - 28 September 2013
c
c################################################################
      integer i,j,k
      
      GCSCavity = tGCSCavity
      if (boundary.eq.PERIODICBOUNDARIES) then
        do i=1,3
          BoxDimension(i) = tBoxDimension(i)
          BoxOrigin(i) = tBoxOrigin(i)
        enddo
      endif
      
c     now write out the coordinates of each protein
      do i=1,NProteins     
        do j=1,NZAtoms(i)
          do k=1,3
            ZCoords(i,j,k) = tZCoords(i,j,k)
          enddo
        enddo
      enddo

c     now write the coordinates of all of the solutes
      do i=1,NSolutes
        do j=1,3
          do k=1,3
            SolDummy(i,j,k) = tSolDummy(i,j,k)
          enddo
          SolGeomCentre(i,j) = tSolGeomCentre(i,j)
        enddo
        do j=1,NSolAtoms(i)
          do k=1,3
            ZSolCoords(i,j,k) = tZSolCoords(i,j,k)
          enddo
        enddo
      enddo

c     now write the coordinates of all of the solvents
      if (NSolvents.gt.0) then
        do i=1,NSolvents
          do j=1,3
            SvnCentre(i,j) = tSvnCentre(i,j)
          enddo
          do j=1,NSvnAtoms(i)
            do k=1,3
              ZSvnCoords(i,j,k) = tZSvnCoords(i,j,k)
            enddo
          enddo
        enddo
      endif
      
c       now write the coordinates of all of the GC-solutes
      if (NGCSolutes.gt.0) then
        do i=1,NGCSolutes
         thetaparam(i) = tthetaparam(i)
          do j=1,3
            GCSCentre(i,j) = tGCSCentre(i,j)
          enddo       
          do j=1,NGCSAtoms(i)
            do k=1,3
              ZGCSCoords(i,j,k) = tZGCSCoords(i,j,k)
            enddo
          enddo
        enddo
      endif

c     now recalculate the total energy
      GotInitialTotalEnergy = .false.
      call checkGotTotal

      end subroutine