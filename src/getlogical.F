      logical function getLogical(word,deflt)
      implicit none
c#################################################
c
c     This function gets the logical variable
c     from the word.
c
c     If there is no variable
c     then the value in 'deflt' is returned
c
c     (C) Christopher Woods
c
c#################################################

      character*(*) word
      character*4 tmp
      logical deflt

      read(word,*,err=9700,end=9700) getLogical
      return 
           
9700  continue
c     could not get a logical, so see if we can
c     interpret the word ourselves
      tmp = word(1:4)
      call lowercase(tmp)
      
      if (tmp(1:3).eq.'yes' .or.
     .    tmp(1:1).eq.'1' .or.
     .    tmp(1:2).eq.'on' .or.
     .    word(1:4).eq.'true') then
      
        getLogical = .true.
      else if (tmp(1:2).eq.'no' .or.
     .         tmp(1:1).eq.'0' .or.
     .         tmp(1:3).eq.'off' .or.
     .         tmp(1:4).eq.'fals') then
        getLogical = .false.
      else
        getLogical = deflt
      endif
      
      return
      end
