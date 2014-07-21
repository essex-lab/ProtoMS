      integer function lineLength(printstring)
      implicit none
c###################################################################
c
c       This function calculates the length of the string 
c       'printstring'. The length is last position of 
c       a text character + 1
c
c       (this leaves an extra space!)
c
c       (C) Christopher Woods        9 August 2002
c
c###################################################################

      character*(*) printstring
      integer lastpos,i,strlen
                                  
      lastpos = 0
      strlen = len(printstring)
      do i=strlen,1,-1
        if (printstring(i:i).ne.' ') then
          lastpos = i
          goto 101
        endif
      enddo

c     the string is all spaces - return 0
      lineLength = 0
      return

101   continue
                           
      lineLength = lastpos
        
      return
      end
