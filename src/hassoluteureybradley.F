      logical function hasSoluteUreyBradley(isol,zat1,zat2,zat3)
      include 'dimensions.inc'
      include 'solutes.inc'
c###########################################################
c
c     This function returns whether there is a Urey-Bradley term
c     between solute atoms with zmatrix id's zat1,zat2 in the solute
c     isol
c
c     (C) Christopher Woods
c
c###########################################################

      integer zat1,zat2,zat3,isol,nubys,i
      
      hasSoluteUreyBradley = .false.
      
      nubys = NSolUreyBradleys(isol)
      if (nubys.le.0) return
      
      hasSoluteUreyBradley = .true.
      do i=1,nubys
        if (SolUreyBradley(isol,i,1).eq.zat1 .and. SolUreyBradley(isol,i,2).eq.zat2
     .                 .and. SolUreyBradley(isol,i,3).eq.zat3) return
        if (SolUreyBradley(isol,i,1).eq.zat3 .and. SolUreyBradley(isol,i,2).eq.zat2
     .                 .and. SolUreyBradley(isol,i,3).eq.zat3) return
      enddo
      hasSoluteUreyBradley = .false.
      return
      end
