      subroutine pushBackSoluteCoords(isol)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'stacks.inc'
      include 'solutes.inc'
c#################################################################
c
c     This routine saves the coordinates of solute isol
c     in the solute coordinate stack (lifo)
c
c     (C) Christopher Woods    28 February 2002
c
c#################################################################

      integer isol,i,j,nats,inbnds,inangs,indihs,inubys

      if (isol.le.0 .or. isol.gt.MAXSOLUTES) return
      
      nats = NSolAtoms(isol)
      if (nats.le.0) return
      
      solCoordDepth = solCoordDepth + 1
      if (solCoordDepth.gt.MAXSOLCOORDDEPTH) then
        write(printstring,*) "Exceeded solute coordinate stack depth! (",
     .                   MAXSOLCOORDDEPTH,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Increase MAXSOLCOORDDEPTH and recompile")
        call closeProgram(-1,"ProtoMS exiting as the stack is too small")
      endif
      
c     now save the coordinates
      do j=1,3
        do i=1,nats
          stackSol(solCoordDepth,i,j) = ZSolCoords(isol,i,j)
          stackSolZ(solCoordDepth,i,j) = ZSolVal(isol,i,j)
        enddo
c       save the centre and principal axes
        stackSolCentre(solCoordDepth,j) = SolCentre(isol,j)
        stackSolGeomCentre(solCoordDepth,j) = SolGeomCentre(isol,j)
        do i=1,3
          stackSolDummies(solCoordDepth,i,j) = SolDummy(isol,i,j)
        enddo
      enddo

      stackSolRadius(solCoordDepth) = SolRadius(isol)

c     save the bond/angle/dihedral sizes
      inbnds = NSolBonds(isol)
      if (inbnds.gt.0) then
        do i=1,inbnds
          stackSolBonds(solCoordDepth,i) = SolBondLength(isol,i)
        enddo
      endif
      inangs = NSolAngles(isol)
      if (inangs.gt.0) then
        do i=1,inangs
          stackSolAngles(solCoordDepth,i) = SolAngleSize(isol,i)
        enddo
      endif
      indihs = NSolDihedrals(isol)
      if (indihs.gt.0) then
        do i=1,indihs
          stackSolDihedrals(solCoordDepth,i) = SolDihedralSize(isol,i)
        enddo
      endif
      inubys = NSolUreyBradleys(isol)
      if (inubys.gt.0) then
        do i=1,inubys
          stackSolUreyBradleys(solCoordDepth,i) = SolUreyBradleySize(isol,i)
        enddo
      endif

c     save the solute number
      stackSolCoord(solCoordDepth) = isol
      return
      end
