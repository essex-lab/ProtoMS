      subroutine UpdateSplitMoveResidueList()
      include 'dimensions.inc'
      include 'proteins.inc'
      include 'quickpot.inc'
      include 'movelist.inc'
      include 'printer.inc'
c################################################
c
c
c
c     
c###############################################
      integer i,iprot,ires
      logical inresiduesplitmoved

c     for every moved residue...
      do i=1,nResMoved
        iprot = ResMovedList(i,1)
        ires = ResMovedList(i,2)
c       check if it is not already in the splitlist
        if ( .not. inresiduesplitmoved(iprot,ires)) then
c         Add it to the SplitMovedList...
c          write(printstring,*) 'Res ',ires,' ',iprot,' must be saved '
c          call printLine(DEBUG,printstring)
          nResSplitMoved = nResSplitMoved + 1
          if (nResSplitMoved .gt. MAXSPLITRESMOVED) then 
            write(printstring,*) "Exceeded split residue list length! (",
     .        MAXSPLITRESMOVED,")"
            call printLine(WARNING,printstring)
            call printLine(WARNING,"Increase MAXSPLITRESMOVED and recompile")
            call closeProgram(-1,"ProtoMC exiting as the list is too small")
          endif
c         Save Residue
          SplitResMovedList(nResSplitMoved,1) = iprot
          SplitResMovedList(nResSplitMoved,2) = ires
        else
c         write(printstring,*) 'Residue ',ires,' ',iprot,
c     .      ' already in split list'
c          call printLine(DEBUG,printstring)
        endif
      enddo
      return
      end
