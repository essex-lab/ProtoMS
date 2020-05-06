      subroutine readInfoLine(adel,rdel)
      include 'dimensions.inc'
      include 'splitter.inc'
c######################################################
c
c     This routine scans the words to find out whether 
c     there are 'rotate' or 'translate' directives.
c
c     (C) Christopher Woods, September 2004
c
c######################################################

      double precision adel,rdel
      double precision readDouble
      double precision deg2Rad
      integer i
      
c     scan through the words to find 'translate' or 'rotate'
      do i=1,nwords
        if (words(i).eq.'translate') then
          adel = readDouble(words(i+1))
        else if (words(i).eq.'rotate') then
          rdel = deg2Rad(readDouble(words(i+1)))
        endif
      enddo
      
      end
