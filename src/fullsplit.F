      subroutine fullsplit(line,char,startpos,maxsplit,maxlen,array,nsplit)
      implicit none
c###################################################################
c
c     This routine will split the line 'line' by character 'char'
c     starting from position 'startpos' in the line. The results
c     are returned in a character array ('array') where each
c     string in the array has maximum length 'maxlen', and the 
c     array holds a maximum of 'maxlen' strings. This routine
c     will also match the number of split strings 'nsplit'.
c
c     e.g. if line = 'mary had a   little lamb', startpos = 3,
c          char = ' ', maxsplit = 6 and maxlen = 4, then
c          the returned array will contain;
c
c          1: 'ry  '  2: 'had '  3: 'a   '  4: 'litt'  5: 'lamb'  6: undefined
c   
c          nsplit would equal 5. Note how multiple matches of the split
c          character are ignored (i.e. between 'a' and 'little'), 
c          how the matched strings are truncated into the limit of 
c          maxlen ('little' is shortened to 'litt'), how the matching
c          started at position 3, so missed of the 'ma' of 'mary', and
c          how only 5 matches were found, so the 6th element of the 
c          returned array is undefined.
c
c     This routine is very complicated and flaky, mainly because
c     fortran should not be used to do this!
c
c     However, it does work! Given
c              1234567890123456789012345      
c      line = "mary had a little lamb, little lamb, little lamb"
c      
c      call split(line," ",3,6,4,array,nsplit)
c
c      returns... Matches are 'ry  ' 'had ' 'a   ' 'litt' 'lamb' 'litt'
c
c     (C) Christopher Woods, 9 September 2004
c
c###################################################################

c     arguments passed to subroutine
      character*(*) line
      character*1 char
      integer startpos
      integer maxsplit
      integer maxlen
      integer nsplit
      character*(*) array(*)

c     subroutine variables
      integer lngth,i
      integer begstr,endstr,endmax
      integer lineLength

      nsplit = 0

c     get the length of the line to be split
      lngth = lineLength(line)

      if (lngth.le.0) return
      
c     make sure that we don't include the last newline character
      if (line(lngth:lngth).eq."\n") then
        lngth = lngth - 1
      endif
           
c     make sure that we are starting at a valid point on the line
      if (startpos.gt.lngth .or. startpos.le.0) return    

c     the first split must start from startpos
      begstr = startpos

c     loop over all of the characters in the line
      do i=startpos,lngth
c       if this character matches the split character...
        if (line(i:i).eq.char .or. i.eq.lngth) then
c         the end of the string is here!
          endstr = i
c         if begstr == endstr then this is an empty match (unless this is a single
c         character match at the end of the string!)
          if (begstr.ne.endstr) then
c           increment the number of matches
            nsplit = nsplit + 1
c           make sure that we don't exceed maxsplit
            if (nsplit.gt.maxsplit) then
               nsplit = maxsplit
               return
            endif
            
c           The matched string will go from begstr to endstr-1.
c           Make sure that the length of this string is less than maxlen
            if (endstr - begstr .gt. maxlen) then
              endmax = begstr + maxlen - 1
            else
              if (line(i:i).eq.char) then
                endmax = endstr - 1
              else 
                endmax = endstr
              endif
            endif
           
c           copy this portion of the string into the array
            array(nsplit) = line(begstr:endmax)
          else if (i.eq.lngth .and. line(i:i).ne.char) then
            nsplit = nsplit + 1
            array(nsplit) = line(i:i)
          endif            
c         now increment the start of the next matched string
          begstr = endstr + 1

        endif
      enddo

c     clear the rest of the array      
      do i=nsplit+1,maxsplit
        array(i) = " "
      enddo
      
      return
      end
