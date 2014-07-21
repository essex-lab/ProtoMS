
      integer function readInt(word)
      include 'dimensions.inc'
      include 'printer.inc'
c##########################################################
c
c     This function reads an integer from
c     a word and returns it
c
c     (C) Christopher Woods, September 2004
c
c##########################################################

      character*(*) word
      
      read(word,*,err=9700,end=9800) readInt
      return 
           
9700  continue
      write(printstring,10) word
10    format("Error trying to read an integer from '",a,"'")
      call printLine(WARNING,printstring)
      readInt = 0
      return
      
9800  continue
      write(printstring,20) word
20    format("Reached end of word '",a,"' while trying to read an integer")
      call printLine(WARNING,printstring)
      readInt = 0
      
      return
      end
