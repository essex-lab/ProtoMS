      subroutine split2(line,token)
      include 'dimensions.inc'
      include 'splitter.inc'
      include 'printer.inc'
c###################################################################
c
c     This routine is similar to split, except split2 will split 
c     by whatever character you want, and return the results in the array
c     words2
c
c     Julien Michel February 2005
c
c###################################################################

      character*(*) line
      character*1 token
      integer i
      
      call fullsplit(line,token,1,SPLITTER_MAXSPLIT,SPLITTER_MAXLEN,
     .                          words2,nwords2)
c     if any of the words start with '#' then truncate the line there and then
      if (nwords2.gt.0) then
        do i=1,nwords2
          if (words2(i)(1:1).eq.'#') then
            nwords2 = i-1
            return
          endif
        enddo
      endif

      return
      end
