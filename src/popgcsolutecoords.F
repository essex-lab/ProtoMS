      subroutine popGCSoluteCoords
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'solvents.inc'
      include 'stacks.inc'
c##############################################################
c
c     This routine pops the coordinates of the last gc-solute
c     from the coordinate stack (lifo)
c
c     MSB 02/04/09
c###############################################################

      integer i,j,nats,igcs
      
      if (gcsCoordDepth.le.0) then
        call printLine(WARNING,"We have popped off the end of the gc-solute coordinate stack!")
        call printLine(WARNING,"This is a serious programming bug!")
        call closeProgram(-1,"ProtoMC is exiting as it has detected a serious programming bug!")
      endif

      igcs = stackGCSCoord(gcsCoordDepth)
      
      if (igcs.le.0 .or. igcs.gt.NGCSolutes) then
        call printLine(WARNING,"The gc-solute stack contains an invalid solute molecule!")
        call printLine(WARNING,"This is a serious programming bug!")
        call closeProgram(-1,"ProtoMC is exiting as it has detected a serious programming bug!")
      endif
      
      nats = NGCSAtoms(igcs)
      if (nats.le.0) then
        call printLine(WARNING,"The solute stack contains an invalid number of atoms!")
        call printLine(WARNING,"This is a serious programming bug!")
        call closeProgram(-1,"ProtoMC is exiting as it has detected a serious programming bug!")
      endif
      
      do j=1,3
        do i=1,nats
          ZGCSCoords(igcs,i,j) = stackGCS(gcsCoordDepth,i,j)
        enddo
        GCSCentre(igcs,j) = stackGCSCentre(gcsCoordDepth,j)
      enddo

      gcsCoordDepth = gcsCoordDepth - 1
      return
      end
