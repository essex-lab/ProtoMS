      subroutine split(line)
      include 'dimensions.inc'
      include 'splitter.inc'
c###################################################################
c
c     This routine is a wrapper around 'fullsplit' to simplify the 
c     splitting of lines where the splitting uses spaces.
c     The results are placed into the splitter common block
c     (include splitter.inc)
c
c     (C) Christopher Woods,  September 2004
c
c###################################################################

      character*(*) line
      integer i
      
      call fullsplit(line," ",1,SPLITTER_MAXSPLIT,SPLITTER_MAXLEN,
     .                          words,nwords)
      
c     if any of the words start with '#' then truncate the line there and then
      if (nwords.gt.0) then
        do i=1,nwords
          if (words(i)(1:1).eq.'#') then
            nwords = i-1
            return
          endif
        enddo
      endif

      return
      end
