      subroutine takeAcid(iprot,ires)
      include 'dimensions.inc'
      include 'proteins.inc'
c################################################################
c
c     This removes an acid from the acid array, and shuffles the
c     rest along to fill the hole...
c
c     (C) Christopher Woods  4 February 2002
c
c################################################################

      integer iprot,ires,i
      logical foundres
      
      foundres = .false.
      
      do i=1,NAcidicRes(iprot)
        if (AcidicRes(iprot,i).eq.ires) foundres = .true.
        if (foundres .and. i.le.MAXACIDICRES) 
     .              AcidicRes(iprot,i) = AcidicRes(iprot,i+1)
      enddo

      if (foundres) then
        NAcidicRes(iprot) = NAcidicRes(iprot) - 1
        if (NAcidicRes(iprot).lt.0) NAcidicRes(iprot) = 0
      endif
      
      return
      end
