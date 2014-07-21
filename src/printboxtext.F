      subroutine printBoxText(STREAM,line)
      include 'dimensions.inc'
      include 'splitter.inc'
      
      integer wdth
      integer i,nspc
      integer STREAM
      integer lineLength,lgth
      integer blnce
      double precision ratio
      integer nwds,istrt,iend
      character*1 lftchar
      character*1 rgtchar
      character*300 boxstr,tmpstr
      character*(*) line
      lftchar = "<"
      rgtchar = ">"
      wdth = 96

      lgth = lineLength(line)
      if (lgth.le.0) then
        write(boxstr,10) lftchar,(" ",i=1,wdth-2),rgtchar
        call printLine(STREAM,boxstr)
      else
c       break the line into individual lines of length < wdth-4
        call split(line)
c       what is the ratio of the width of the box and the width of the line?
        ratio = dble(wdth-4) / dble(lgth)
c       number of words per line is ratio * nwords
        if (ratio.ge.1.D+00) then
          nwds = nwords
        else
          nwds = int(ratio*nwords)
        endif
        
c       now loop until all of the words are written        
        istrt = 0
        iend = 0
100     continue
          istrt = iend + 1
          iend = istrt + nwds - 1

          if (istrt.gt.nwords) return
          if (iend.gt.nwords) iend = nwords
          
101       continue          
          call join(istrt,iend,tmpstr)

          lgth = lineLength(tmpstr)
          if (lgth.gt.wdth-4) then
            iend = iend - 1
            if (iend.eq.istrt) then
c             the word is longer than width! just print all of it out on this line!              
              iend = iend + 1
            else
              goto 101
            endif
          endif
               
          nspc = ( (wdth-lgth)/2 ) - 1
          if (nspc.lt.1) nspc = 1
          
          blnce = wdth-lgth-2-(2*nspc)
          if (blnce.lt.0) blnce = 0
          
          write(boxstr,10) lftchar,(" ",i=1,nspc),tmpstr(1:lgth),(" ",i=1,nspc+blnce),rgtchar
10        format(300a)
          call printLine(STREAM,boxstr)
        goto 100
      endif
      
      return
      end
