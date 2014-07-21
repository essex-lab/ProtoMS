
      subroutine clearSurfaces()
      include 'dimensions.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'printer.inc'
      include 'simparams.inc'
      include 'surface.inc'
c######################################
c     Reset the Old and New surfaces of 
c     the entire loaded system
c
c     Julien Michel November 2004
c#######################################
      integer i,j,k,nres,nats,resstart,resend

      if (NProteins .gt. 0) then
        do i=1,NProteins
          nres = NResidues(i)
          do j=1,nres
            OldResidueSurf(i,j) = ZERO
            OldResidueSurfF(i,j) = ZERO
            OldResidueSurfB(i,j) = ZERO
            NewResidueSurf(i,j) = ZERO
            NewResidueSurfF(i,j) = ZERO
            NewResidueSurfB(i,j) = ZERO
            resstart = BBZid(i,j,1)
            resend = BBZid(i,j,4)
            do k=resstart,resend
              OldProtAtomSurf(i,k) = ZERO
              OldProtAtomSurfF(i,k) = ZERO
              OldProtAtomSurfB(i,k) = ZERO
              NewProtAtomSurf(i,k) = ZERO
              NewProtAtomSurfF(i,k) = ZERO
              NewProtAtomSurfB(i,k) = ZERO
            enddo
          enddo
        enddo
      endif
      
      if (NSolutes .gt. 0) then
        do i=1,NSolutes
          OldSoluteSurf(i) = ZERO
          OldSoluteSurfF(i) = ZERO
          OldSoluteSurfB(i) = ZERO
          NewSoluteSurf(i) = ZERO
          NewSoluteSurfF(i) = ZERO
          NewSoluteSurfB(i) = ZERO
          nats = NSolAtoms(i)
          do j=1,nats
            OldSolAtomSurf(i,j) = ZERO
            OldSolAtomSurfF(i,j) = ZERO
            OldSolAtomSurfB(i,j) = ZERO
            NewSolAtomSurf(i,j) = ZERO
            NewSolAtomSurfF(i,j) = ZERO
            NewSolAtomSurfB(i,j) = ZERO
          enddo
        enddo
      endif

      if (NSolvents .gt. 0) then
        do i=1,NSolvents
          OldSolventSurf(i) = ZERO
          OldSolventSurfF(i) = ZERO
          OldSolventSurfB(i) = ZERO
          NewSolventSurf(i) = ZERO
          NewSolventSurfF(i) = ZERO
          NewSolventSurfB(i) = ZERO
          nats = NSvnAtoms(i)
          do j=1,nats
            OldSvnAtomSurf(i,j) = ZERO
            OldSvnAtomSurfF(i,j) = ZERO
            OldSvnAtomSurfB(i,j) = ZERO
            NewSvnAtomSurf(i,j) = ZERO
            NewSvnAtomSurfF(i,j) = ZERO
            NewSvnAtomSurfB(i,j) = ZERO
          enddo
        enddo
      endif

      return
      end
