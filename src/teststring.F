      logical function testString(string)
c#################################################
c     this function tests to make sure 
c     that the string contains input
c     (is not empty)
c#################################################
      character*(*) string
      integer lineLength
      integer length

      length = lineLength(string)
      
      if (length.gt.0) then
        teststring = .true.
      else
        teststring = .false.
      endif
      
      end
