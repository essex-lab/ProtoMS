      subroutine SaveSolventGB(isvn)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'printer.inc'
      include 'gb.inc'
c############################################
c     Save solvent coordinates on to a list
c     of moving solvents for GB calculations
c
c     Julien Michel November 2004
c##############################################
      integer isvn,nats,i,j

      if (isvn.le.0 .or. isvn.gt.MAXSOLVENTS) return
      nats = NSvnAtoms(isvn)
      if (nats.le.0) return
      
      countGBSolvent = countGBSolvent + 1
      if (countGBSolvent .gt. MAXSOLVENTS) then
        call printLine(WARNING,"Increase MAXSOLVENTS  and recompile")
        call closeProgram(-1,"ProtoMC exiting !")
      endif
      
c     now save the coordinates
      do j=1,3
        do i=1,nats
          gbSolventCoords(countGBSolvent,i,j) = ZSvnCoords(isvn,i,j)
        enddo
c       save the centre and principal axes
        gbSvnCentre(countGBSolvent,j) = SvnCentre(isvn,j)
      enddo
c     Save solute number
      gbSolventList(countGBSolvent) = isvn
      
      return
      end


