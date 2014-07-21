      subroutine SaveSoluteGB(isol)
      include 'dimensions.inc'
      include 'solutes.inc'
      include 'printer.inc'
      include 'gb.inc'
c############################################
c     Save solute coordinates on to a list
c     of moving solutes for GB calculations
c
c     Julien Michel November 2004
c##############################################
      integer isol,nats,i,j

      if (isol.le.0 .or. isol.gt.MAXSOLUTES) return
      nats = NSolAtoms(isol)
      if (nats.le.0) return
      
      countGBSolute = countGBSolute + 1
      if (countGBSolute .gt. MAXSOLUTES) then
        call printLine(WARNING,"Increase MAXSOLUTES  and recompile")
        call closeProgram(-1,"ProtoMC exiting !")
      endif
      
c     now save the coordinates
      do j=1,3
        do i=1,nats
          gbSoluteCoords(countGBSolute,i,j) = ZSolCoords(isol,i,j)
          if ( NSolVariable(isol) .gt. 0 ) then
            gbSoluteCoordsF(countGBSolute,i,j) = ZSolCoords(NSolForwards(isol),i,j)
            gbSoluteCoordsB(countGBSolute,i,j) = ZSolCoords(NSolBackwards(isol),i,j)
          endif
        enddo
c       save the centre and principal axes
        gbSolCentre(countGBSolute,j) = SolCentre(isol,j)
      enddo
      gbSolRadius(countGBsolute) = SolRadius(isol)
c     Save solute number
      gbSoluteList(countGBSolute) = isol
      
      return
      end


