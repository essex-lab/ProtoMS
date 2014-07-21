      subroutine join(istrt,iend,line)
      include 'dimensions.inc'
      include 'splitter.inc'
c#########################################################
c
c     Join together words istrt to iend and return
c     the result in 'line'
c
c     (C) Christopher Woods, October 2004
c
c#########################################################

      integer istrt,iend,idel,i
      character*(*) line
      integer lgth,lineLength,strtpos,endpos,l
      
      if (istrt.le.0) istrt = 1
      if (iend.le.0) iend = 1
      if (istrt.gt.nwords) istrt = nwords
      if (iend.gt.nwords) iend = nwords
      
      if (istrt.gt.iend) then
        idel = -1
      else
        idel = 1
      endif
      
      strtpos = 0
      endpos = 0
      lgth = len(line)
      do i=istrt,iend,idel
        l = lineLength(words(i))
        endpos = strtpos + l
        if (endpos.gt.lgth) then
          endpos = lgth
          l = endpos - strtpos
          if (l.eq.0) return
        endif
      
        strtpos = strtpos + 1
        line(strtpos:endpos) = words(i)(1:l)
        strtpos = endpos + 1
        if (strtpos.gt.lgth) return
        line(strtpos:strtpos) = " "
      enddo
      
      do i=strtpos+1,lgth
        line(i:i) = " "
      enddo
      
      return
      end
