      subroutine pushBackSolventCoords(isvn)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'printer.inc'
      include 'enums.inc'
      include 'stacks.inc'
c########################################################
c
c     This routine saves the coordinates of solvent
c     isvn in the solvent coordinate stack
c
c     (C) Christopher Woods   5 March 2002
c
c########################################################

      integer isvn,i,j,nats

      if (isvn.le.0 .or. isvn.gt.NSolvents) return
      
      nats = NSvnAtoms(isvn)
      if (nats.le.0) return
      
      svnCoordDepth = svnCoordDepth + 1
      if (svnCoordDepth.gt.MAXSVNCOORDDEPTH) then
        write(printstring,*) "Exceeded solvent coordinate stack depth! (",
     .                   MAXSVNCOORDDEPTH,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Increase MAXSVNCOORDDEPTH and recompile")
        call closeProgram(-1,"ProtoMC is exiting as the solvent stack is too small")
      endif
      
c     now save the coordinates, as well as the centre
      do j=1,3
        do i=1,nats
          stackSvn(svnCoordDepth,i,j) = ZSvnCoords(isvn,i,j)
        enddo
        StackSvnCentre(svnCoordDepth,j) = SvnCentre(isvn,j)
c        StackSvnAxis(SvnCoordDepth,1,j) = SvnMajAxis(isvn,j)
c        StackSvnAxis(SvnCoordDepth,2,j) = SvnMinAxis(isvn,j)
c        StackSvnAxis(SvnCoordDepth,3,j) = SvnPerAxis(isvn,j)
      enddo

c     save the solvent number
      stackSvnCoord(svnCoordDepth) = isvn
      return
      end
