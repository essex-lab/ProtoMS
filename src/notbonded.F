      logical function notBonded(ityp,ival,zat1,zat2,bond14)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'solutes.inc'
      include 'proteins.inc'
c#############################################################
c
c     This function returns whether or not atoms with
c     zid's zat1 and zat2 in ityp/ival are not 1-1, 1-2 or 1-3 bonded. If
c     they are 1-4 bonded, bond14 is true
c
c     (C) Christopher Woods
c
c#############################################################

      integer ityp,ival,iprot,isol,i,ires,ires2,nbnds,nangs,ndihs,zat1,zat2
      integer ibnd,iang,idih
      logical bond14
      
      bond14 = .false.
      notBonded = .false.
      
      if (zat1.eq.zat2) return
      
      if (ityp.eq.PROTEIN) then
        iprot = ival
        ires = Zid2ResID(iprot,zat1,1)
        ires2 = Zid2ResID(iprot,zat2,1)
      
c       check to see if the residues are bonded - if not, then they cannot
c       have any valid 1-4 interactions...
        if (ires.ne.ires2) then
          if (.not.ResResBond(iprot,ires,ires2)) then
            bond14 = .true.
            notBonded = .true.
            return
          endif
        endif
     
        nbnds = NResBonds(iprot,ires)
        if (nbnds.gt.0) then
          do i=1,nbnds
            ibnd = ResidueBond(iprot,ires,i)
            if (Bond(iprot,ibnd,1).eq.zat1 .and. Bond(iprot,ibnd,2).eq.zat2) return
            if (Bond(iprot,ibnd,1).eq.zat2 .and. Bond(iprot,ibnd,2).eq.zat1) return
          enddo
        endif
       
        nangs = NResAngles(iprot,ires)
        if (nangs.gt.0) then
          do i=1,nangs
            iang = ResidueAngle(iprot,ires,i)
            if (Angle(iprot,iang,1).eq.zat1 .and. Angle(iprot,iang,3).eq.zat2) return
            if (Angle(iprot,iang,1).eq.zat2 .and. Angle(iprot,iang,3).eq.zat1) return
          enddo
        endif
      
        notBonded = .true.
        ndihs = NResDihedrals(iprot,ires)
        if (ndihs.gt.0) then
          do i=1,ndihs
            idih = ResidueDihedral(iprot,ires,i)
            if ((Dihedral(iprot,idih,1).eq.zat1 .and. Dihedral(iprot,idih,4).eq.zat2) .or.
     .        (Dihedral(iprot,idih,1).eq.zat2 .and. Dihedral(iprot,idih,4).eq.zat1)) then
               bond14 = .true.
               return
            endif
          enddo
        endif
      
        return

      else if (ityp.eq.SOLUTE) then
        isol = ival
        nbnds = NSolBonds(isol)
        if (nbnds.gt.0) then
          do i=1,nbnds
            if (SolBond(isol,i,1).eq.zat1 .and. SolBond(isol,i,2).eq.zat2) return
            if (SolBond(isol,i,1).eq.zat2 .and. SolBond(isol,i,2).eq.zat1) return
          enddo
        endif
        nangs = NSolAngles(isol)
        if (nangs.gt.0) then
          do i=1,nangs
              if (SolAngle(isol,i,1).eq.zat1 .and. SolAngle(isol,i,3).eq.zat2) return
              if (SolAngle(isol,i,1).eq.zat2 .and. SolAngle(isol,i,3).eq.zat1) return
          enddo
        endif
        ndihs = NSolDihedrals(isol)
        notBonded = .true.
        if (ndihs.gt.0) then
          do i=1,ndihs
            if ((SolDihedral(isol,i,1).eq.zat1 .and. SolDihedral(isol,i,4).eq.zat2)
     .      .or.(SolDihedral(isol,i,1).eq.zat2 .and. SolDihedral(isol,i,4).eq.zat1)) then
              bond14 = .true.
              return
            endif
          enddo
        endif
        
        return
        
      endif        
      
      end
