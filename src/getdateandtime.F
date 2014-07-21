      subroutine getDateAndTime(year,month,day,hour,mins,secs)
      implicit none
c#####################################################################
c
c     This function gets the current date and time and returns
c     it in the obviously named integers.
c
c     (C) Christopher Woods, October 2004
c
c#####################################################################

      integer year,month,day,hour,mins,secs
      character*15 datestr,timestr
      
      year = 1970
      month = 1
      day = 1
      hour = 0
      mins = 0
      secs = 0
      
      call Date_and_Time(datestr,timestr)

c     datestr has the format yyyymmdd
      read(datestr,10) year,month,day
10    format(i4,i2,i2)      
            
c     timestr has the format hhmmss.xx (xx are milliseconds)
      read(timestr,11) hour,mins,secs
11    format(i2,i2,i2)      

      return
      end
