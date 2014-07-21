      subroutine UpdateSplitMoveSolventList()
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'quickpot.inc'
      include 'movelist.inc'
      include 'printer.inc'
c################################################
c
c
c
c     
c###############################################
      integer i,isvn
      logical insolventsplitmoved

c     for every moved solute...
      do i=1,nSvnMoved
        isvn = SvnMovedList(i)
c       check if it is not already in the splitlist
        if ( .not. insolventsplitmoved(isvn)) then
c         Add it to the SplitMovedList...
c          write(printstring,*) 'Solvent ',isvn,' must be saved '
c          call printLine(DEBUG,printstring)
          nSvnSplitMoved = nSvnSplitMoved + 1
          if (nSvnSplitMoved .gt. MAXSPLITSVNMOVED) then 
            write(printstring,*) "Exceeded split solute list length! (",
     .        MAXSPLITSVNMOVED,")"
            call printLine(WARNING,printstring)
            call printLine(WARNING,"Increase MAXSPLITSOLMOVED and recompile")
            call closeProgram(-1,"ProtoMC exiting as the list is too small")
          endif
c         Save Solute
          SplitSvnMovedList(nSvnSplitMoved) = isvn
c        else
c          write(printstring,*) 'Solvent ',isvn,' already in split list'
c          call printLine(DEBUG,printstring)
        endif
      enddo
      return
      end
