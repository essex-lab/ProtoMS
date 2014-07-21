      logical function hasUreyBradley(iprot,zat1,zat2,zat3)
      include 'dimensions.inc'
      include 'proteins.inc'
c#################################################
c
c     This function returns true is protein iprot
c     contains a Urey-Bradley between z atoms zat1-zat2-zat3
c
c#################################################

      integer iprot,zat1,zat2,zat3
      integer ires,nubys,iuby,i
      
      hasUreyBradley = .false.
      
      ires = Zid2ResID(iprot,zat1,1)
      
      nubys = NResUreyBradleys(iprot,ires)
      if (nubys.le.0) return
      
      hasUreyBradley = .true.
      do i=1,nubys
        iuby = ResidueUreyBradley(iprot,ires,i)
        if (UreyBradley(iprot,iuby,1).eq.zat1 .and. 
     .      UreyBradley(iprot,iuby,2).eq.zat2 .and.
     .      UreyBradley(iprot,iuby,3).eq.zat3) return
        if (UreyBradley(iprot,iuby,1).eq.zat3 
     .      .and. UreyBradley(iprot,iuby,2).eq.zat2
     .      .and. UreyBradley(iprot,iuby,3).eq.zat1) return
      enddo

      hasUreyBradley = .false.
      return
      end
