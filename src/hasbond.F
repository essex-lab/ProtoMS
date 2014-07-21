      logical function hasBond(iprot,zat1,zat2)
      include 'dimensions.inc'
      include 'proteins.inc'
c#################################################
c
c     This function returns true is protein iprot
c     contains a bond between z atoms zat1 and zat2
c
c#################################################

      integer iprot,zat1,zat2
      integer ires,nbnds,ibnd,i
      
      hasBond = .false.
      
      ires = Zid2ResID(iprot,zat1,1)
      
      nbnds = NResBonds(iprot,ires)
      if (nbnds.le.0) return
      
      hasBond = .true.
      do i=1,nbnds
        ibnd = ResidueBond(iprot,ires,i)
        if (Bond(iprot,ibnd,1).eq.zat1 .and. Bond(iprot,ibnd,2).eq.zat2) return
        if (Bond(iprot,ibnd,1).eq.zat2 .and. Bond(iprot,ibnd,2).eq.zat1) return
      enddo

      hasBond = .false.
      return
      end
