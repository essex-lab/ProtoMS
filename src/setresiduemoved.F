      subroutine setResidueMoved(iprot,ires)
      include 'dimensions.inc'
      include 'proteins.inc'
      include 'movelist.inc'
      include 'printer.inc'
c#################################################
c
c     This sets the residue 'ires' in protein
c     'iprot' as having moved. Note that this
c     does *not* check to see if this residue
c     is already in the list!
c
c     (C) Christopher Woods, October 2004
c
c#################################################

      integer iprot,ires,in
      
c     do not need to do this if everything is already moved
      if (allMoved) return

      if (iprot.le.0 .or. iprot.gt.NProteins) return
      if (ires.le.0 .or. ires.gt.NResidues(iprot)) return
      
      in = nResMoved + 1
      if (in.gt.MAXRESMOVED) then
        write(printstring,*) "Exceeded maximum number of moved residues (",MAXRESMOVED,")"
        call printLine(FATAL,printstring)
        write(printstring,*)"Increase MAXRESMOVED and recompile"
        call closeProgram(-1,printstring)
        return
      endif
      
      resMovedList(in,1) = iprot
      resMovedList(in,2) = ires
      
      nResMoved = in
      return
      end
