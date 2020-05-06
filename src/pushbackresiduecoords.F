      subroutine pushBackResidueCoords(iprot,ires)
      include 'dimensions.inc'
      include 'stacks.inc'
      include 'proteins.inc'
      include 'printer.inc'
c#######################################################
c
c     This routine saves the coordinates of residue ires 
c     in protein iprot in the residue coordinates stack
c
c     (C) Christopher Woods   11 February 2002
c
c########################################################

      integer i,j,iprot,ires,nsc
      
      if (iprot.le.0 .or. iprot.gt.NProteins) return
      if (ires.le.0 .or. ires.gt.NResidues(iprot)) return
      
      coordDepth = coordDepth + 1
      if (coordDepth.gt.MAXRESCOORDDEPTH) then 
        write(printstring,*) "Exceeded residue coordinate stack depth! (",
     .       MAXRESCOORDDEPTH,") - increase MAXRESCOORDDEPTH and recompile"
        call closeProgram(-1,printstring)
      endif
      
c     now save the old residue coordinates
      nsc = NSCAtoms(iprot,ires)
      do i=1,3
        do j=1,4
          stackResBB(coordDepth,j,i) = BBCoord(iprot,ires,j,i)
        enddo
         
        if (nsc.gt.0) then
          do j=1,nsc
            stackResSC(coordDepth,j,i) = SCCoord(iprot,ires,j,i)
            stackResZSC(coordDepth,j,i) = ZValSC(iprot,ires,j,i)
          enddo
        endif
        stackResCentre(coordDepth,i) = ResCentre(iprot,ires,i)
      enddo

      stackResRadius(coordDepth) = ResRadius(iprot,ires)
      stackResCoord(coordDepth) = ires
      stackResProt(coordDepth) = iprot

      return
      end
