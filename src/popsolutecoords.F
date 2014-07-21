      subroutine popSoluteCoords
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'solutes.inc'
      include 'stacks.inc'
c############################################################
c
c     This routine pops the coordinates of a solute off
c     from the solute coordinate stack (lifo)
c
c     (C) Christopher Woods   28 February 2002
c
c############################################################

      integer i,j,nats,isol,inbnds,inangs,indihs,inubys
      
      if (solCoordDepth.le.0) then
        call printLine(WARNING,"We have popped off the end of the solute coordinate stack!")
        call printLine(WARNING,"This is a serious programming bug!")
        call closeProgram(-1,"ProtoMC is exiting as it has detected a serious programming bug!")
      endif
      
      isol = stackSolCoord(solCoordDepth)
      
      if (isol.le.0 .or. isol.gt.NSolutes) then
        call printLine(WARNING,"The solute stack contains an invalid solute!")
        call printLine(WARNING,"This is a serious programming bug!")
        call closeProgram(-1,"ProtoMC is exiting as it has detected a serious programming bug!")
      endif
      
      nats = NSolAtoms(isol)
      if (nats.le.0) then
        call printLine(WARNING,"The solute stack contains an invalid number of atoms!")
        call printLine(WARNING,"This is a serious programming bug!")
        call closeProgram(-1,"ProtoMC is exiting as it has detected a serious programming bug!")
      endif
      
      do j=1,3
        do i=1,nats
          ZSolCoords(isol,i,j) = stackSol(solCoordDepth,i,j)
          ZSolVal(isol,i,j) = stackSolZ(solCoordDepth,i,j)
        enddo
        SolCentre(isol,j) = stackSolCentre(solCoordDepth,j)
        SolGeomCentre(isol,j) = stackSolGeomCentre(solCoordDepth,j)
        do i=1,3
          SolDummy(isol,i,j) = stackSolDummies(solCoordDepth,i,j)
        enddo
      enddo

      SolRadius(isol) = stackSolRadius(solCoordDepth)

c     restore the bond/angle/dihedral sizes      
      inbnds = NSolBonds(isol)
      if (inbnds.gt.0) then
        do i=1,inbnds
          SolBondLength(isol,i) = stackSolBonds(solCoordDepth,i)
        enddo
      endif
      inangs = NSolAngles(isol)
      if (inangs.gt.0) then
        do i=1,inangs
          SolAngleSize(isol,i) = stackSolAngles(solCoordDepth,i)
        enddo
      endif
      indihs = NSolDihedrals(isol)
      if (indihs.gt.0) then
        do i=1,indihs
          SolDihedralSize(isol,i) = stackSolDihedrals(solCoordDepth,i)
        enddo
      endif
      inubys = NSolUreyBradleys(isol)
      if (inubys.gt.0) then
        do i=1,inubys
          SolUreyBradleySize(isol,i) = stackSolUreyBradleys(solCoordDepth,i)
        enddo
      endif
      
      solCoordDepth = solCoordDepth - 1
      return
      end
