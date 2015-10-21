      subroutine shutdown
      use timing
      include 'printer.inc'
      include 'version.inc'
c###########################################################
c
c     This routine shuts down the simulation and prints
c     out the results and restart information
c
c     (C) Christopher Woods     9 August 2002
c
c###########################################################

c     variables that hold the current date and time
      integer year,month,day,hour,mins,secs

      write(printstring,*) "Simulation shutting down..."
      call printLine(INFO,printstring)

c     print the time that everything finished
      call getDateAndTime(year,month,day,hour,mins,secs)
      write(printstring,10) exename,hour,mins,secs,day,month,year
10    format(a," completed at ",i2,":",i2.2,":",i2.2," on the ",i2,"/",i2,"/",i4)      
      call printLine(INFO,printstring)

      call timer_end(EMAIN)
      call print_timing
            
      return
      end
