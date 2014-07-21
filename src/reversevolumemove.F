      subroutine reverseVolumeMove
      include 'dimensions.inc'
      include 'proteins.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'solvents.inc'
      include 'solutes.inc'
      include 'pressure.inc'
      include 'enums.inc'
      include 'flags.inc'
c##########################################################
c
c       This subroutine reverses the previous volume
c       move - it can only be called directly after
c       a volume move, or weird things will happen!
c
c       (C) Christopher Woods  30 April 2003
c
c##########################################################

      integer i,j,k,l
      double precision idel

      if (boundary .eq. PERIODICBOUNDARIES) then
        do i=1,3
          BoxDimension(i) = OldBoxDimension(i)
          BoxOrigin(i) = OldBoxOrigin(i)
          BoxTop(i) = BoxOrigin(i) + BoxDimension(i)
          BoxDimension2(i) = 0.5D+00*BoxDimension(i)
        enddo
        volume = BoxDimension(1)*BoxDimension(2)*BoxDimension(3)          

        if (NProteins.gt.0) then
          do i=1,NProteins
            do k=1,3
              idel = OldProCentre(i,k)-ProCentre(i,k)
              if (idel.ne.ZERO) then
                do j=1,NZAtoms(i)
                  ZCoords(i,j,k) = ZCoords(i,j,k) + idel
                enddo
                ProCentre(i,k) = ProCentre(i,k) + idel
                do j=1,NResidues(i)
                  ResCentre(i,j,k) = ResCentre(i,j,k)+idel
                  do l=1,4
                    BBCoord(i,j,l,k) = BBCoord(i,j,l,k) + idel
                  enddo
                  do l=1,NSCAtoms(i,j)
                    SCCoord(i,j,l,k) = SCCoord(i,j,l,k) + idel
                  enddo
                enddo
              endif
            enddo
          enddo
        endif
        
        if (NSolutes.gt.0) then
          do i=1,NSolutes
            do k=1,3
              idel = OldSolCentre(i,k)-SolCentre(i,k)
              if (idel.ne.ZERO) then
                do j=1,NSolAtoms(i)
                  ZSolCoords(i,j,k) = ZSolCoords(i,j,k) + idel
                enddo
                SolCentre(i,k) = OldSolCentre(i,k)
                SolGeomCentre(i,k) = OldSolGeomCentre(i,k)
              endif
            enddo
          enddo
        endif

      else if (boundary .eq. SPHERICALHARMONIC) then
        capradius = oldcapradius
      endif

      if (NSolvents.gt.0) then
        do i=1,NSolvents
          do k=1,3
            idel = OldSvnCentre(i,k)-SvnCentre(i,k)
            if (idel.ne.ZERO) then
              do j=1,NSvnAtoms(i)
                ZSvnCoords(i,j,k) = ZSvnCoords(i,j,k) + idel
              enddo
              SvnCentre(i,k) = OldSvnCentre(i,k)
            endif
          enddo
        enddo
      endif
        
      return
      end
