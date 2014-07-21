      subroutine assignSurface
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'constants.inc'
      include 'printer.inc'
      include 'surface.inc'
c#################################
c     Assign the loaded Surface Area terms 
c     to each CLJ type that was loaded
c
c##################################
      integer i,j,findSA

      write(printstring,*) 'Assign Surface Area Terms'
      call printLine(INFO,printstring)

      do 10 i=1,MAXCLJ
        if (AtmDB(i) .eq. NULLCLJ) goto 10
        j = findSA(AtmDB(i))
        if (j .lt. 0) then
          write(printstring,*) ' Could not assign radius/surface '
     &      ," tension to AMBER atom ",AtmDB(i)," setting radius "
     &      ," to ZERO. This is probably not what you want !"
          call printLine(WARNING,printstring)
          SurfRadDB(i) = 0.0d0
          SurfTensionDB(i) = 0.0d0
        else
          SurfRadDB(i) = SARadius(j)
          SurfTensionDB(i) = SASurftens(j)          
c          write(printstring,*) ' ipar ',i,' at ',AtmDB(i),
c     &      ' now has a surface ',SurfRadDB(i),' tens ',SurfTensionDB(i)
c          call printLine(DEBUG,printstring)
        endif
 10   continue

      return
      end

