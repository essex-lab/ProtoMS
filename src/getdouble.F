      double precision function getDouble(word,deflt)
      implicit none
c#################################################
c
c     This function gets the real value
c     from the word.
c
c     If there is no variable
c     then the value in 'deflt' is returned
c
c     (C) Christopher Woods
c
c#################################################

      character*(*) word
      double precision deflt
      
      read(word,*,err=9700,end=9700) getDouble
      
      return 
           
9700  continue
c     could not get a double, so set default value
      getDouble = deflt
      
      return

      end
