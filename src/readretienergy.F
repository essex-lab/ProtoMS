      subroutine readRETIenergy(filename,lam,lamf,lamb,nrg,nrgf,nrgb,stat)
      include 'dimensions.inc'
      include 'reti.inc'
      include 'splitter.inc'
      include 'printer.inc'
c###################################################################
c
c     This routine reads the energy from a RETI energy file
c     written by the RETIEnergy routine
c
c     Samuel Genheden, 2012
c
c###################################################################    
      integer FILE 
      parameter (FILE=87)
      integer stat,lineLength
      character*255 filename,line
      double precision lam,lamf,lamb,nrg,nrgf,nrgb,getDouble
      
      stat = 0
      open(unit=FILE,file=filename,status="old",form='formatted',err=9600)

      read(FILE,'(A255)')line
      call split(line)
      if (nwords.ne.7) then
        write(printstring,"(5A,I2,A)")"Corrupt RETI file(",filename(1:lineLength(filename)),") : ",words(1)," ",nwords," ! Exit."
        call prettyLine(FATAL,printstring) 
        stat = 9500
        return   
      endif
c      write(printstring,*)"Line to split: ",line," ",nwords
c      call prettyLine(INFO,printstring)

      lam = getDouble(words(2),0.0d0)
      lamb = getDouble(words(4),0.0d0)
      lamf = getDouble(words(6),0.0d0)
      
      nrg = getDouble(words(3),0.0d0)
      nrgb = getDouble(words(5),0.0d0)
      nrgf = getDouble(words(7),0.0d0)

      if ((lam.eq.0.0d0).and.(lamb.eq.0.0d0).and.(lamf.eq.0.0d0)) then
        write(printstring,"(3A)")"Non-sense lambda values in RETI file (",filename(1:lineLength(filename)),") ! Exit."
        call prettyLine(FATAL,printstring) 
        stat = 9400
        return   
      endif

      return
      
 9600 write(printstring,"(3A)")"Unable to open RETI file (",filename(1:lineLength(filename)),") ! Exit."
      call prettyLine(FATAL,printstring) 
      stat = 9600    
      return 

      end subroutine

