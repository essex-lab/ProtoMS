      subroutine revertCoordinates(conf)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'solvents.inc'
      include 'flags.inc'
      include 'energies.inc'
      include 'printer.inc'
      include 'quickpot.inc'      
c########################################
c
c     This subroutine copies back in place
c     the coordinates of a selected 
c     configuration
c
c     (c) Julien Michel 2005
c########################################
      integer i,j,k,conf

c      write(printstring,*) 'REVERT COORDINATES OF CONF ',conf
c      call printLine(DEBUG,printstring)

      if (NProteins .gt. 0) then
        do i=1,NProteins
          do j=1,NZAtoms(i)
            do k=1,3
              ZCoords(i,j,k) = StoreZCoords(i,j,k,conf) 
            enddo
          enddo
        enddo
      endif

      if (NSolutes.gt.0) then
        do i=1,NSolutes
          do k=1,3
            do j=1,NSolAtoms(i)
              ZSolCoords(i,j,k) = StoreZSolCoords(i,j,k,conf) 
            enddo
            SolCentre(i,k) = StoreSolCentre(i,k,conf) 
            SolGeomCentre(i,k) = StoreSolGeomCentre(i,k,conf) 
            do j=1,3
              SolDummy(i,j,k) = StoreSolDummies(i,j,k,conf) 
            enddo
          enddo
        enddo
      endif

      if (NSolvents .gt. 0) then
        do i=1,NSolvents
          do k=1,3
            do j=1,NSvnAtoms(i)
              ZSvnCoords(i,j,k) = StoreZSvnCoords(i,j,k,conf) 
            enddo
            SvnCentre(i,k) = StoreSvnCentre(i,k,conf)
          enddo
        enddo
      endif

c
c     Now rebuild internals
c
      if (NSolutes.gt.0) then
        do i=1,NSolutes
          call getSoluteZmatrix(i)
          call updateSoluteInternals(i)
          call getSoluteCentreradius(i)
        enddo
      endif

      if (NProteins .gt. 0) then
        do i=1,NProteins
          call syncProteinZmatrix(i)
        enddo
      endif


      return
      end
