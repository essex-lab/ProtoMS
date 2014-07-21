      subroutine pushBackGCSoluteCoords(igcs)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'printer.inc'
      include 'enums.inc'
      include 'stacks.inc'
c########################################################
c
c     This routine saves the coordinates of gcsolute
c     igcs in the gc-solute coordinate stack
c
c     (C) Christopher Woods   5 March 2002
c
c########################################################

      integer igcs,i,j,nats

      if (igcs.le.0 .or. igcs.gt.NGCSolutes) return
      
      nats = NGCSAtoms(igcs)
      if (nats.le.0) return
      
      gcsCoordDepth = gcsCoordDepth + 1
      if (gcsCoordDepth.gt.MAXSVNCOORDDEPTH) then
        write(printstring,*) "Exceeded solvent coordinate stack depth! (",
     .                   MAXSVNCOORDDEPTH,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Increase MAXSVNCOORDDEPTH and recompile")
        call closeProgram(-1,"ProtoMC is exiting as the solvent stack is too small. Test msb")
      endif
      
c     now save the coordinates, as well as the centre
      do j=1,3
        do i=1,nats
          stackGCS(gcsCoordDepth,i,j) = ZGCSCoords(igcs,i,j)
        enddo
        StackGCSCentre(gcsCoordDepth,j) = GCSCentre(igcs,j)
c        StackSvnAxis(SvnCoordDepth,1,j) = SvnMajAxis(isvn,j)
c        StackSvnAxis(SvnCoordDepth,2,j) = SvnMinAxis(isvn,j)
c        StackSvnAxis(SvnCoordDepth,3,j) = SvnPerAxis(isvn,j)
      enddo

c     save the solvent number
      stackGCSCoord(gcsCoordDepth) = igcs
      return
      end
