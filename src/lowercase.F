      subroutine lowercase(line)
      implicit none
c###################################################################
c
c     This function converts the contents of the string 'line' to
c     all lowercase.
c
c     (C) Christopher Woods, September 2004
c
c###################################################################

      character*(*) line
      integer lgth,i
      
      lgth = len(line)

      do i=1,lgth
        if (line(i:i).ge.'A' .and. line(i:i).le.'Z') then
          line(i:i) = char(ichar(line(i:i))+32)
        endif
      enddo
      
      return
      end
