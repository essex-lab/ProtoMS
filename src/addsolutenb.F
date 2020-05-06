      subroutine addSoluteNB(isol,i,j,is14)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'solutes.inc'
c#############################################################
c
c     Quick subroutine to add a non-bonded interaction
c     to solute isol between atoms i and j, and with
c     is14=true if it is a 14 interaction
c
c     (C) Christopher Woods
c
c#############################################################

      integer i,j,isol
      logical is14
      
      NSolNB(isol) = NSolNB(isol) + 1

      if (NSolNB(isol).gt.MAXSOLUTENBPAIRS) then
        write(printstring,*) "Exceeded maximum number of solute nb pairs (",
     .                        MAXSOLUTENBPAIRS,")"
        call printLine(FATAL,printstring)
        write(printstring,*)"Increase MAXSOLUTENBPAIRS and recompile"
        call closeProgram(-1,printstring)
        NSolNB(isol) = MAXSOLUTENBPAIRS
        return
      endif

      SolNB(isol,NSolNB(isol),1) = i
      SolNB(isol,NSolNB(isol),2) = j
      SolNB14(isol,NSolNB(isol)) = is14

      return
      end
