      logical function hasSoluteBond(isol,zat1,zat2)
      include 'dimensions.inc'
      include 'solutes.inc'
c###########################################################
c
c     This function returns whether there is a bond between
c     solute atoms with zmatrix id's zat1,zat2 in the solute
c     isol
c
c     (C) Christopher Woods   19 Febrary 2002
c
c###########################################################

      integer zat1,zat2,isol,nbnds,i
      
      hasSoluteBond = .false.
      
      nbnds = NSolBonds(isol)
      if (nbnds.le.0) return
      
      hasSoluteBond = .true.
      do i=1,nbnds
        if (SolBond(isol,i,1).eq.zat1 .and. SolBond(isol,i,2).eq.zat2) return
        if (SolBond(isol,i,1).eq.zat2 .and. SolBond(isol,i,2).eq.zat1) return
      enddo
      hasSoluteBond = .false.
      return
      end
