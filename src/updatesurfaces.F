      
      subroutine updateSurfaces()
      include 'dimensions.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'solvents.inc'
      include 'energies.inc'
      include 'flags.inc'
      include 'constants.inc'
      include 'surface.inc'
      include 'printer.inc'
c############################################
c
c
c     Julien Michel November 2004
c############################################
      integer i,j,k,nats,nres,resstart,resend
        
      if (NSolutes .gt. 0) then
        do i=1,NSolutes
          if (SoluteredoSA(i)) then
            OldSoluteSurf(i) = NewSoluteSurf(i) 
            OldSoluteSurfF(i) = NewSoluteSurfF(i)
            OldSoluteSurfB(i) = NewSoluteSurfB(i)
            nats = NSolAtoms(i)
            do j=1,nats
              if (SolAtredoSA(i,j)) then
                OldSolAtomSurf(i,j) = NewSolAtomSurf(i,j)
                OldSolAtomSurfF(i,j) = NewSolAtomSurfF(i,j)
                OldSolAtomSurfB(i,j) = NewSolAtomSurfB(i,j)
              endif
            enddo
          endif
        enddo
      endif

      if (NProteins .gt. 0) then
        do i=1,NProteins
          OldProteinSurf(i) = NewProteinSurf(i)
          OldProteinSurfF(i) = NewProteinSurfF(i)
          OldProteinSurfB(i) = NewProteinSurfB(i)
          nres = NResidues(i)
          do j=1,nres
            if (ResidueredoSA(i,j)) then
              OldResidueSurf(i,j) = NewResidueSurf(i,j)
              OldResidueSurfF(i,j) = NewResidueSurfF(i,j)
              OldResidueSurfB(i,j) = NewResidueSurfB(i,j)
              resstart = BBZid(i,j,1)
              resend = BBZid(i,j,4)
              do k=resstart,resend
                if (ProtAtredoSA(i,k)) then
                  OldProtAtomSurf(i,k) = NewProtAtomSurf(i,k)
                  OldProtAtomSurfF(i,k) = NewProtAtomSurfF(i,k)
                  OldProtAtomSurfB(i,k) = NewProtAtomSurfB(i,k)
                endif
              enddo
            endif
          enddo
        enddo
      endif

      if (NSolvents .gt. ZERO) then
        do i=1,NSolvents
          if (SolventredoSA(i)) then
            OldSolventSurf(i) = NewSolventSurf(i)
            OldSolventSurfF(i) = NewSolventSurfF(i)
            OldSolventSurfB(i) = NewSolventSurfB(i)
            nats = NSvnAtoms(i)
            do j=1,nats
              if (SvnAtredoSA(i,j)) then
                OldSvnAtomSurf(i,j) = NewSvnAtomSurf(i,j)
                OldSvnAtomSurfF(i,j) = NewSvnAtomSurfF(i,j)
                OldSvnAtomSurfB(i,j) = NewSvnAtomSurfB(i,j)
              endif
            enddo
          endif
        enddo
      endif
      
      return
      end
