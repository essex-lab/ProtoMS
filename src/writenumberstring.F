      subroutine writeNumberString(str,val,reslt)
      implicit none
c#################################################
c
c     This function attaches the integer number 'val'
c     to the end of the string 'str', e.g.
c
c     writeNumberString('protein',32) will return
c     a string that is 'protein32' in reslt
c
c     (C) Christopher Woods
c
c#################################################

       character*(*) str
       character*(*) reslt
       integer val,q,r
       integer i,npow
       integer diff,strlen
       integer lineLength

c      use a character array to map individual integers to
c      characters       
       character*1 ints(10)
       data ints /'0','1','2','3','4','5','6','7','8','9'/
       
       strlen = lineLength(str)
       
c      make sure we don't overflow the string
       diff = len(reslt) - strlen
       if (diff.le.0 .or. val.lt.0) then
         if (strlen.gt.len(reslt)) strlen = len(reslt)
         reslt = str(1:strlen)
         return
       endif

c      copy the string
       reslt(1:strlen) = str(1:strlen)

       npow = 0
c      find the highest power of ten in the string - only go up to 10^diff
       do i=1,diff
         q = val / 10**i
         if (q.eq.0) goto 100
         npow = i
       enddo
       
       reslt = str(1:len(reslt))
       return
       
c      get here only if the number can fit into the string
100    continue
       
c      now divide the string by 10^pow and get the quotients
       r = val
       do i=npow,1,-1
         q = r / 10**i
         r = r - q*10**i
         strlen = strlen + 1
         reslt(strlen:strlen) = ints(q+1)
       enddo

c      now add on the last remainder
       strlen = strlen + 1
       reslt(strlen:strlen) = ints(r+1)
       strlen = strlen + 1

c      now fill the rest of the string with spaces
       if (strlen.lt.len(reslt)) then
         do i=strlen,len(reslt)
           reslt(i:i) = ' '
         enddo
       endif       
       
       return
       end
