      subroutine assignBorn
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'constants.inc'
      include 'printer.inc'
      include 'gb.inc'
c#################################
c     Assign the loaded Generalised Born terms 
c     to each CLJ type that was loaded
c
c##################################
      integer i,j,findGB

      write(printstring,*) 'Assigning Generalised Born Terms'
      call printLine(INFO,printstring)

      do 10 i=1,MAXCLJ
        if (AtmDB(i) .eq. NULLCLJ) goto 10
        j = findGB(AtmDB(i))
        if (j .lt. 0) then
          write(printstring,*) ' Could not assign born parameters'
     &      ," to AMBER atom ",AtmDB(i)," setting them   "
     &      ," to ZERO. This is probably not what you want !"
          call printLine(WARNING,printstring)
          IntrBornDB(i) = 0.0d0
          PDAScaleDB(i) = 0.0d0
        else
          IntrBornDB(i) = MAX(SigDB(i) * SIGTORS * GBIBorn(j),MINIBORN)
          PDAScaleDB(i) = IntrBornDB(i) * GBScale(j)
c          write(printstring,*) ' ipar ',i,' at ',AtmDB(i),
c     &      ' now has a  ',IntrBornDB(i),' ',PDAScaleDB(i)
c          call printLine(DEBUG,printstring)
        endif
 10   continue

      return
      end

