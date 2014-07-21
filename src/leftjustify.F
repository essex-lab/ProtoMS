      subroutine leftJustify(line)
      implicit none
c################################################
c
c     This routine left justifies the contents
c     of line
c
c################################################

      character*(*) line
      integer lgth,lineLength,i,strt
      
      lgth = lineLength(line)
      if (lgth.le.0) return
      
      strt = 0
      
c     find the start of the words
      do i=1,lgth
        if (line(i:i).ne.' ') then
          strt = i
          goto 100
        endif
      enddo
100   continue

c     if strt is 0 or 1 then nothing needs to be done
      if (strt.le.1) return
      
      do i=strt,lgth
        line(i-strt+1:i-strt+1) = line(i:i)
      enddo
      
      do i=lgth-strt+2,lgth
        line(i:i) = ' '
      enddo
      
      return
      end
