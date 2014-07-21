      logical function hasSoluteDihedral(isol,zat1,zat2,zat3,zat4)
      include 'dimensions.inc'
      include 'solutes.inc'
c############################################################
c
c     This function returns whether there is a dihedral angle between
c     solute atoms zat1-zat2-zat3-zat4 in solute isol
c
c     (C) Christopher Woods   19 February 2002
c
c############################################################
      
      integer isol,zat1,zat2,zat3,zat4,ndihs,i

      hasSoluteDihedral = .false.
      ndihs = NSolDihedrals(isol)
      if (ndihs.le.0) return
      hasSoluteDihedral = .true.
      do i=1,ndihs
        if (SolDihedral(isol,i,1).eq.zat1 .and. SolDihedral(isol,i,2).eq.zat2
     .  .and. SolDihedral(isol,i,3).eq.zat3 .and. SolDihedral(isol,i,4).eq.zat4
     .   ) return
        if (SolDihedral(isol,i,1).eq.zat4 .and. SolDihedral(isol,i,2).eq.zat3
     .  .and. SolDihedral(isol,i,3).eq.zat2 .and. SolDihedral(isol,i,4).eq.zat1
     .   ) return
      enddo
      hasSoluteDihedral = .false.
      return
      end
