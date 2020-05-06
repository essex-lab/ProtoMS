      subroutine initialise
      include 'dimensions.inc'
      include 'printer.inc'
      include 'version.inc'
c###########################################################
c
c     This routine initialises the program
c
c     (C) Christopher Woods     9 August 2002
c
c###########################################################

c     variables that hold the current date and time
      integer year,month,day,hour,mins,secs
      
C     print the program header
      call printHeader

c     print the time that everything started      
      call getDateAndTime(year,month,day,hour,mins,secs)
      write(printstring,10) exename,hour,mins,secs,day,month,year
10    format(a," started at ",i2,":",i2.2,":",i2.2," on the ",i2,"/",i2,"/",i4)      
      call printLine(INFO,printstring)
      
C     call routine to initialise all the data
      call initialiseArrays

      end
