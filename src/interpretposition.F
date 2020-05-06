
      integer function interpretPosition(word)
      implicit none
      include 'enums.inc'
c###################################################
c
c     This function interprets the word to 
c     work out which protein position it refers
c     to. It then returns the protein position
c
c     (C) Christopher Woods, November 2004
c
c###################################################

      character*(*) word

      if (word.eq.'first') then
        interpretPosition = FIRST
      else if (word.eq.'middle') then
        interpretPosition = MIDDLE
      else if (word.eq.'last') then
        interpretPosition = LAST
      else if (word.eq.'single') then
        interpretPosition = SINGLE
      else if (word.eq.'all') then
        interpretPosition = ALLPOS
      else
        interpretPosition = NOPOS
      endif
      
      return
      end
      
