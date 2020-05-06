      subroutine printBoxBar(STREAM)
      implicit none
      
      integer wdth
      integer i
      integer STREAM
      character*1 lftchar
      character*1 rgtchar
      character*200 boxstr
      
      lftchar = "<"
      rgtchar = ">"
      wdth = 96
      
      write(boxstr,10) (lftchar,rgtchar,i=1,wdth/2)
10    format(200A)

      call printLine(STREAM,boxstr)
      
      return
      end
