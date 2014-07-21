      subroutine takeBase(iprot,ires)
      include 'dimensions.inc'
      include 'proteins.inc'
c################################################################
c
c     This removes an base from the base array, and shuffles the
c     rest along to fill the hole...
c
c     (C) Christopher Woods  4 February 2002
c
c################################################################

      integer iprot,ires,i
      logical foundres
      
      foundres = .false.
      
      do i=1,NBasicRes(iprot)
        if (BasicRes(iprot,i).eq.ires) foundres = .true.
        if (foundres .and. i.le.MAXBASICRES) 
     .              BasicRes(iprot,i) = BasicRes(iprot,i+1)
      enddo

      if (foundres) then
        NBasicRes(iprot) = NBasicRes(iprot) - 1
        if (NBasicRes(iprot).lt.0) NBasicRes(iprot) = 0
      endif
      
      return
      end
