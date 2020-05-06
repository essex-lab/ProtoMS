      subroutine getDOFInfo(dummy,flex,delta,ipar0,ipar1)
      include 'dimensions.inc'
      include 'splitter.inc'
c##################################################################
c
c     This routine scans the words to find out whether a 
c     degree of freedom is a dummy, whether it is flexible,
c     and if so what the maximum move delta is, and whether
c     there are any overriding parameters (-1 if there aren't)
c
c     (C) Christopher Woods, September 2004
c
c##################################################################

      logical dummy,flex
      double precision delta
      integer ipar0,ipar1
      integer i
      integer readInt
      double precision readDouble
      
      flex = .false.
      dummy = .false.
      delta = 0.D+00
      ipar0 = -1
      ipar1 = -1
      
c     scan through the words to find 'dummy', 'flex' or 'param'      
      do i=1,nwords
        if (words(i).eq.'dummy') then
          dummy = .true.
        else if (words(i).eq.'flex') then
c         format is 'flex idelta'
          flex = .true.
          delta = readDouble(words(i+1))
        else if (words(i).eq.'param') then
c         format is 'param ipar0 ipar1'
          ipar0 = readInt(words(i+1))
          ipar1 = readInt(words(i+2))
        
c         make sure that parameters less than 0 are interpreted as dummies
          if (ipar0.lt.0) ipar0 = 0
          if (ipar1.lt.0) ipar1 = 0
        endif
      enddo
      
      return
      end
