      subroutine StoreCoordinates(conf)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'solvents.inc'
      include 'flags.inc'
      include 'energies.inc'
      include 'printer.inc'
      include 'quickpot.inc'
c############################################################
c     This subroutine save the coordinates of the entire 
c     system
c
c     (c) Julien Michel April 2005
c###########################################################
      integer i,j,k,conf

c      write(printstring,*) 'STORING COORDINATES of conf ',conf
c      call printLine(DEBUG,printstring)

      if (NProteins .gt. 0) then
        do i=1,NProteins
          do j=1,NZAtoms(i)
            do k=1,3
              StoreZCoords(i,j,k,conf) = ZCoords(i,j,k)
            enddo
          enddo
        enddo
      endif

      if (NSolutes.gt.0) then
        do i=1,NSolutes
          do k=1,3
            do j=1,NSolAtoms(i)
              StoreZSolCoords(i,j,k,conf) = ZSolCoords(i,j,k)
            enddo
            StoreSolCentre(i,k,conf) = SolCentre(i,k)
            StoreSolGeomCentre(i,k,conf) = SolGeomCentre(i,k)
            do j=1,3
              StoreSolDummies(i,j,k,conf) = SolDummy(i,j,k)
            enddo
          enddo
        enddo
      endif

      if (NSolvents .gt. 0) then
        do i=1,NSolvents
          do k=1,3
            do j=1,NSvnAtoms(i)
              StoreZSvnCoords(i,j,k,conf) = ZSvnCoords(i,j,k)
            enddo
            StoreSvnCentre(i,k,conf) = SvnCentre(i,k)
          enddo
        enddo
      endif

      return
      end
