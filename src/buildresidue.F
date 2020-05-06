      subroutine buildResidue(iprot,ires)
      include 'dimensions.inc'
      include 'proteins.inc'
c##################################################
c
c     This subroutine rebuilds residue ires from 
c     iprot, from the backbone atom coordinates,
c     and the current values in the zmatrix...
c
c     It also ensures that the coordinates in both
c     the residue based (BBatom and SCatom) and zmatrix-indexed
c     based (ZCoords) arrays are in synch.
c
c     (C) Christopher Woods, October 2004
c
c##################################################

      integer iprot,ires
      integer i,j,natoms

c     go over all the sidechain atoms and build those first...
      natoms = NSCAtoms(iprot,ires)
      if (natoms.gt.0) then
        do i=1,natoms
          call genProteinAtom(iprot,ZMatSC(iprot,ires,i,3),ZMatSC(iprot,ires,i,2),
     .                              ZMatSC(iprot,ires,i,1),SCZid(iprot,ires,i),
     .                              ZValSC(iprot,ires,i,3),ZValSC(iprot,ires,i,2),
     .                              ZValSC(iprot,ires,i,1))
          
c         duplicate these coordinates from the zmatrix-array to the residue-array
          do j=1,3
             SCCoord(iprot,ires,i,j) = ZCoords(iprot,SCZid(iprot,ires,i),j)
          enddo
        enddo
      endif
      
c     duplicate the atoms to the residue-based coordinate array     
      do i=1,4
        do j=1,3
          BBCoord(iprot,ires,i,j) = ZCoords(iprot,BBZid(iprot,ires,i),j)
        enddo
      enddo
c     update the residue internals to reflect the new coordinates
      call updateResidueInternals(iprot,ires)
c     get the new residue centre and radius
      call getResidueCentreRadius(iprot,ires)
      return 
     
      end
