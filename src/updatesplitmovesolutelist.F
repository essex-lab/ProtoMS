      subroutine UpdateSplitMoveSoluteList()
      include 'dimensions.inc'
      include 'solutes.inc'
      include 'quickpot.inc'
      include 'movelist.inc'
      include 'printer.inc'
c################################################
c
c
c
c     
c###############################################
      integer i,isol
      logical insolutesplitmoved

c     for every moved solute...
      do i=1,nSolMoved
        isol = SolMovedList(i)
c       check if it is not already in the splitlist
        if ( .not. insolutesplitmoved(isol)) then
c         Add it to the SplitMovedList...
c          write(printstring,*) 'Solute ',isol,' must be saved '
c          call printLine(DEBUG,printstring)
          nSolSplitMoved = nSolSplitMoved + 1
          if (nSolSplitMoved .gt. MAXSPLITSOLMOVED) then 
            write(printstring,*) "Exceeded split solute list length! (",
     .        MAXSPLITSOLMOVED,")"
            call printLine(WARNING,printstring)
            call printLine(WARNING,"Increase MAXSPLITSOLMOVED and recompile")
            call closeProgram(-1,"ProtoMS exiting as the list is too small")
          endif
c         Save Solute
          SplitSolMovedList(nSolSplitMoved) = isol
c        else
c          write(printstring,*) 'Solute ',isol,' already in split list'
c          call printLine(DEBUG,printstring)
        endif
      enddo
      return
      end
