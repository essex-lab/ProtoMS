
      subroutine moveVolume(newvol)
      include 'dimensions.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'simparams.inc'
      include 'solvents.inc'
      include 'solutes.inc'
      include 'pressure.inc'
      include 'enums.inc'
      include 'flags.inc'
c###################################################
c
c     This routine changes the volume of the
c     system to newvol
c
c     (C) Christopher Woods   30 April 2003
c
c###################################################

      double precision newvol
      double precision icent(3)
      double precision iscl,idel
      double precision imin,imax
      integer i,j,k,l

c     the scale factor is equal to the cube root of the change in volume
      iscl = (newvol/volume)**(ONE/THREE)

      if (boundary .eq. PERIODICBOUNDARIES) then
c       The old dimensions of the solvent box are in 'BoxDimension(3)'
c       (and 0.5*BoxDimension(3) in BoxDimension2(3))

c       need to scale the coordinates *from the centre of the box!*
        do i=1,3
          icent(i) = BoxOrigin(i) + BoxDimension2(i)
          OldBoxOrigin(i) = BoxOrigin(i)
          BoxOrigin(i) = icent(i) + iscl*(BoxOrigin(i)-icent(i))
          BoxTop(i) = icent(i) + iscl*(BoxTop(i)-icent(i))
          OldBoxDimension(i) = BoxDimension(i)
          BoxDimension(i) = iscl*BoxDimension(i)
          BoxDimension2(i) = 0.5*BoxDimension(i)
        enddo

      else if (boundary .eq. SPHERICALHARMONIC) then
c       The old cap radius is in capradius - we need to 
c       scale the radius by 'iscl' (as r_new = r_old * cbrt(v_new/v_old))
        oldcapradius = capradius
        capradius = iscl * capradius

c       need to scale the coordinates *from the centre of the cap!*
        do i=1,3
          icent(i) = cap(i)
        enddo
      endif

c     copy the centres of the molecules to old storage - this will allow the 
c     volume move to be reversed if it is rejected...
      if (NProteins.gt.0) then
        do i=1,NProteins
          do k=1,3
            imin = 999.0
            imax = -999.0
c           need to calculate the centre of the protein - this is the centre of the 
c           residue centres!
            do j=1,NResidues(i)
              if (ResCentre(i,j,k).lt.imin) imin = ResCentre(i,j,k)
              if (ResCentre(i,j,k).gt.imax) imax = ResCentre(i,j,k)
            enddo

            ProCentre(i,k) = imin + 0.5*(imax-imin)
            OldProCentre(i,k) = ProCentre(i,k)
          enddo
        enddo
      endif
        
      if (NSolutes.gt.0) then
        do i=1,NSolutes
          do j=1,3
            OldSolCentre(i,j) = SolCentre(i,j)
            OldSolGeomCentre(i,j) = SolGeomCentre(i,j)
          enddo
        enddo
      endif
        
      if (NSolvents.gt.0) then
        do i=1,NSolvents
          do j=1,3
            OldSvnCentre(i,j) = SvnCentre(i,j)
          enddo
        enddo
      endif
        
c     now we can move everything (well, we only move the 
c     proteins and solutes if we are using periodic boundaries)
      if (boundary .eq. PERIODICBOUNDARIES) then
        if (NProteins.gt.0) then
          do i=1,NProteins
            do k=1,3
c             get the delta for this dimension - scaled from box centre!
              idel = (iscl-ONE)*(ProCentre(i,k)-icent(k))
c             add this delta to all of the coordinates
              ProCentre(i,k) = ProCentre(i,k) + idel
              do j=1,NZAtoms(i)
                ZCoords(i,j,k) = ZCoords(i,j,k) + idel
              enddo
              do j=1,NResidues(i)
                ResCentre(i,j,k) = ResCentre(i,j,k) + idel
                do l=1,4
                  BBCoord(i,j,l,k) = BBCoord(i,j,l,k) + idel
                enddo
                do l=1,NSCAtoms(i,j)
                  SCCoord(i,j,l,k) = SCCoord(i,j,l,k) + idel
                enddo
              enddo
            enddo
          enddo
        endif
             
        if (NSolutes.gt.0) then
          do i=1,NSolutes
            do k=1,3
              idel = (iscl-ONE)*(SolCentre(i,k)-icent(k))
              SolCentre(i,k) = SolCentre(i,k) + idel
              SolGeomCentre(i,k) = SolGeomCentre(i,k) + idel
              do j=1,NSolAtoms(i)
                ZSolCoords(i,j,k) = ZSolCoords(i,j,k) + idel
              enddo
            enddo
          enddo
        endif
      endif
 
      if (NSolvents.gt.0) then
        do i=1,NSolvents
          do k=1,3
            idel = (iscl-ONE)*(SvnCentre(i,k)-icent(k))
            SvnCentre(i,k) = SvnCentre(i,k) + idel
            do j=1,NSvnAtoms(i)
              ZSvnCoords(i,j,k) = ZSvnCoords(i,j,k) + idel
            enddo
          enddo
        enddo
      endif
        
      volume = newvol
        
      return
        
      end
