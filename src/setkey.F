      subroutine setKey(key,value)
      implicit none
      include 'printer.inc'
      include 'keyset.inc'

      character*(*) key,value

c-------------SUBROUTINE setKey(Key,value)-----------------      
C     Called by ParseFile to add key/value pairs to the
c     dictionary
     
      nkeys = nkeys + 1
      
      if (nkeys.gt.MAXKEYS) then
        write(printstring,*) "Exceeded maximum number of keys (",MAXKEYS,")"
        call printLine(FATAL,printstring)
        write(printstring,*)"Increase MAXKEYS and recompile"
        call closeProgram(-1,printstring)
      endif
      
      keys(nkeys) = key
      values(nkeys) = value
      
      return     
      
      end     
