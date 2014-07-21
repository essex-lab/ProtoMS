      integer function getInt(word,deflt)
      implicit none
c#################################################
c
c     This function gets an integer value
c     from the word
c
c     If there is no variable
c     then the value in 'deflt' is returned
c
c     (C) Christopher Woods
c
c#################################################

      character*(*) word
      integer deflt

      read(word,*,err=9700,end=9700) getInt
      return 
           
9700  continue
c     could not get an int, so set default value
      getInt = deflt
      
      return

      end
