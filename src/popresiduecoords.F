      subroutine popResidueCoords
      include 'dimensions.inc'
      include 'stacks.inc'
      include 'proteins.inc'
      include 'printer.inc'
c#####################################################
c
c     This subroutine pops the coordinates of the 
c     residue on the top of the stack (LIFO) back
c     into the residue and zmatrix coordinate arrays
c
c     (C) Christopher Woods    11 February 2002
c
c#####################################################

      integer iprot,ires,i,j,nsc

      if (coordDepth.le.0) then
        write(printstring,*) "We have popped off the end of the coordinate stack! - ",
     .                       "This is a severe programming bug!"
        call closeProgram(-1,printstring)
      endif
      
c     extract the protein and residue number
      iprot = stackResProt(coordDepth)
      ires = stackResCoord(coordDepth)
      
      if (iprot.le.0 .or. iprot.gt.NProteins) then
        call printLine(WARNING,"Serious problem with the coordinate stack!")
        call printLine(WARNING,"This is a result of a programming bug!")
        call closeProgram(-1,"Exiting due to stack bug")
      else if (ires.le.0 .or. ires.gt.NResidues(iprot)) then
        call printLine(WARNING,"Very serious problem with the coordinate stack!")
        call printLine(WARNING,"This is a result of a programming bug!")
        call closeProgram(-1,"Exiting due to stack bug")
      endif
      
c     restore the coordinates to both the residue and zmatrix
c     ordered coordinate arrays
      do i=1,3
        do j=1,4
          BBCoord(iprot,ires,j,i) = stackResBB(coordDepth,j,i)
          ZCoords(iprot,BBZid(iprot,ires,j),i) = BBCoord(iprot,ires,j,i)
        enddo
        
        nsc = NSCAtoms(iprot,ires)
        if (nsc.gt.0) then
          do j=1,nsc
            SCCoord(iprot,ires,j,i) = stackResSC(coordDepth,j,i)
            ZCoords(iprot,SCZid(iprot,ires,j),i) = SCCoord(iprot,ires,j,i)
            ZValSC(iprot,ires,j,i) = stackResZSC(coordDepth,j,i)
          enddo
        endif
        ResCentre(iprot,ires,i) = stackResCentre(coordDepth,i)
      enddo
      
      ResRadius(iprot,ires) = stackResRadius(coordDepth)

      coordDepth = coordDepth - 1

      return
      end
