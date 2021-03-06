      double precision function readDouble(word)
      include 'dimensions.inc'
      include 'printer.inc'
c##########################################################
c
c     This function reads a double precision number from
c     a word and returns it
c
c     (C) Christopher Woods, September 2004
c
c##########################################################

      character*(*) word
      
      read(word,*,err=9700,end=9800) readDouble
      return 
           
9700  continue
      write(printstring,10) word
10    format("Error trying to read a double from '",a,"'")
      call printLine(WARNING,printstring)
      readDouble = 0.0
      return
      
9800  continue
      write(printstring,20) word
20    format("Reached end of word '",a,"' while trying to read a double")
      call printLine(WARNING,printstring)
      readDouble = 0.0
      
      return
      end
