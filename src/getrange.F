      subroutine getRange(word,imin,imax,istrt,iend)
      implicit none
c#########################################################
c
c     This routine tries to read 'word' to work out
c     if it contains a number, or a range of numbers.
c 
c     imin and imax are the min and max possible
c     values of the range. istrt and iend are the
c     returned value of the range. If istrt and iend
c     are both zero, then none of the range is 
c     selected. If imin or imax are less than zero
c     then the range could not be deciphered.
c
c     (C) Christopher Woods, October 2004
c
c########################################################

      character*(*) word
      character*4 tmp
      integer lgth,lineLength,getInt
      integer i,imin,imax,istrt,iend
      
      istrt = -1
      iend = -1
      
      lgth = lineLength(word)
      if (lgth.le.0) return
      
      if (lgth.ge.3) then
        tmp = word(1:4)
        call lowercase(tmp)
        if (tmp(1:3).eq.'all') then
          istrt = imin
          iend = imax
          return
        else if (tmp(1:4).eq.'none') then
          istrt = 0
          iend = 0
          return
        endif
      endif
      
c     now scan each letter for a '-' sign
      if (lgth.ge.3) then
        do i=2,lgth-1
          if (word(i:i).eq.'-') then
c           split the word into two
            istrt = getInt(word(1:i-1),-1)
            iend = getInt(word(i+1:lgth),-1)
            goto 100
          endif
        enddo
      endif

c     we didn't find a '-', so there is only one number            
      istrt = getInt(word,-1)
      iend = istrt
      
100   continue

c     ensure that the range is valid
      if (istrt.lt.0 .or. iend.lt.0) then
        istrt = -1
        iend = -1
        return
      endif
      
      if (istrt.lt.imin) istrt = imin
      if (istrt.gt.imax) istrt = imax
      if (iend.lt.imin) iend = imin
      if (iend.gt.imax) iend = imax
      
      if (iend.lt.istrt) iend = istrt
      
      return
      end
