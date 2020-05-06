      subroutine popSolventCoords
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'solvents.inc'
      include 'stacks.inc'
c##############################################################
c
c     This routine pops the coordinates of the last solvent
c     from the coordinate stack (lifo)
c
c     (C) Christopher Woods   5 March 2002
c
c###############################################################

      integer i,j,nats,isvn
      
      if (svnCoordDepth.le.0) then
        call printLine(WARNING,"We have popped off the end of the solvent coordinate stack!")
        call printLine(WARNING,"This is a serious programming bug!")
        call closeProgram(-1,"ProtoMS is exiting as it has detected a serious programming bug!")
      endif

      isvn = stackSvnCoord(svnCoordDepth)
      
      if (isvn.le.0 .or. isvn.gt.NSolvents) then
        call printLine(WARNING,"The solvent stack contains an invalid solvent molecule!")
        call printLine(WARNING,"This is a serious programming bug!")
        call closeProgram(-1,"ProtoMS is exiting as it has detected a serious programming bug!")
      endif
      
      nats = NSvnAtoms(isvn)
      if (nats.le.0) then
        call printLine(WARNING,"The solvent stack contains an invalid number of atoms!")
        call printLine(WARNING,"This is a serious programming bug!")
        call closeProgram(-1,"ProtoMS is exiting as it has detected a serious programming bug!")
      endif
      
      do j=1,3
        do i=1,nats
          ZSvnCoords(isvn,i,j) = stackSvn(svnCoordDepth,i,j)
        enddo
        SvnCentre(isvn,j) = stackSvnCentre(svnCoordDepth,j)
c        SvnMajAxis(isvn,j) = StackSvnAxis(svnCoordDepth,1,j)
c        SvnMinAxis(isvn,j) = StackSvnAxis(svnCoordDepth,2,j)
c        SvnPerAxis(isvn,j) = StackSvnAxis(svnCoordDepth,3,j)
      enddo

      svnCoordDepth = svnCoordDepth - 1
      return
      end
