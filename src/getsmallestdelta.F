      subroutine getSmallestDelta(istrt,iend)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'constants.inc'
c###############################################################
c
c     This function fixes istrt and iend (two angles in degrees)
c     such that the perturbation between them proceeds via the 
c     smallest possible angle.
c
c     For example, 358 and 2 degrees are converted into 358 and 362
c                  or -2 and 2. Thus the perturbation occurs over
c                  4 degrees, and not 356!
c
c     (C) Christopher Woods   25 January 2002
c
c###############################################################

      double precision istrt,iend
      double precision istrt1,iend1,idel1
      
c     first make sure that both angles lie between 0 and 360
      istrt1 = istrt
      iend1 = iend

      if (istrt1.gt.THREESIXTY) istrt1 = istrt1 - THREESIXTY
      if (istrt1.lt.ZERO) istrt1 = istrt1 + THREESIXTY
      
      if (iend1.gt.THREESIXTY) iend1 = iend1 - THREESIXTY
      if (iend1.lt.ZERO) iend1 = iend1 + THREESIXTY
      
c     now calculate the magnitude of the change
      idel1 = istrt1 - iend1

c     if this is greater than 180 degrees, then there is a smaller
c     perturbation possible
      if (abs(idel1).gt.ONEEIGHTY) then
c       subtract 360 from the angle which is greater than 180              
        if (istrt1.gt.ONEEIGHTY) then
          istrt1 = istrt1 - THREESIXTY
        else if (iend1.gt.ONEEIGHTY) then
          iend1 = iend1 - THREESIXTY
        else
          write(printstring,*) "Something strange in getSmallestDelta(",
     .                     istrt,",",iend,")"
          call printLine(WARNING,printstring)
        endif
        
        idel1 = istrt1 - iend1
        if (abs(idel1).gt.ONEEIGHTY) then
          write(printstring,*) "Something very strange in getSmallestDelta(",
     .          istrt,",",iend,") - ",istrt1," ",iend1
          call printLine(WARNING,printstring)
          return
        endif
        
c       modify the values of istrt and iend accordingly
        istrt = istrt1
        iend = iend1
        return
        
      else
        return
      endif
      
      return
      end
