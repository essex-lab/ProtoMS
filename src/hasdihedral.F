      logical function hasDihedral(iprot,zat1,zat2,zat3,zat4)
      include 'dimensions.inc'
      include 'proteins.inc'
c########################################################
c
c     This function returns whether protein iprot contains
c     a dihedral between zat1-zat2-zat3-zat4
c
c########################################################

      integer iprot,zat1,zat2,zat3,zat4
      integer ndihs,idih,ires,i

      hasDihedral = .false.
      
      ires = Zid2ResID(iprot,zat1,1)
      ndihs = NResDihedrals(iprot,ires)
      if (ndihs.le.0) return
      
      hasDihedral = .true.
      do i=1,ndihs
        idih = ResidueDihedral(iprot,ires,i)
        if (Dihedral(iprot,idih,2).eq.zat2 .and. Dihedral(iprot,idih,3).eq.zat3) then
          if (Dihedral(iprot,idih,1).eq.zat1 .and. Dihedral(iprot,idih,4).eq.zat4) return
        else if (Dihedral(iprot,idih,2).eq.zat3 .and. Dihedral(iprot,idih,3).eq.zat2) then
          if (Dihedral(iprot,idih,1).eq.zat4 .and. Dihedral(iprot,idih,4).eq.zat1) return
        endif
      enddo

      hasDihedral = .false.
      return
      end
