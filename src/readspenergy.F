      subroutine readSPenergy(filename,nrg,stat)
      include 'reti.inc'
      include 'splitter.inc'
      include 'printer.inc'
c###################################################################
c
c     This routine reads the energy from a single point energy file
c     written by the singlePointEnergy routine
c
c     Samuel Genheden, 2013
c
c###################################################################    
      integer FILE 
      parameter (FILE=87)
      integer stat
      character*255 filename,line
      double precision nrg,getDouble
      
      ntries = 0
      stat = 0
      open(unit=FILE,file=filename,status="old",form='formatted',err=9600)

      read(FILE,'(A255)')line
      call split(line)
      if (nwords.ne.8) then
        write(printstring,"(5A,I2,A)")"Corrupt SP file (",filename(1:lineLength(filename)),") : ",words(1)," ",nwords," ! Exit."
        call prettyLine(FATAL,printstring) 
        stat = 9500
        return   
      endif
c      write(printstring,*)"Line to split: ",line," ",nwords
c      call prettyLine(INFO,printstring)

      nrg = getDouble(words(4),0.0d0)

      return
      
 9600 write(printstring,"(3A)")"Unable to open SP file (",filename(1:lineLength(filename)),") ! Exit."
      call prettyLine(FATAL,printstring) 
      stat = 9600    
      return 

      end subroutine

