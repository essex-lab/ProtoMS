      subroutine translateAcid(iprot,ires,vec,del)
      include 'dimensions.inc'
      include 'proteins.inc'
c####################################################
c
c     This function translates the backbone
c     atoms of residue ires in protein iprot
c     by del along the vector vec
c
c     (C) Christopher Woods, October 2004
c
c####################################################
      
      double precision vec(3)
      double precision del,ivec(3)
      integer i,j,iprot,ires,ires1,ires2

c     backbone atoms form a triangle with next residue...
c     (residue A, bb atoms 1-4, residue B, bb atoms 1-4
c      
c              A4          B4
c              |           |
c      A1--A2--A3--B1--B2--B3
c
c      So A2,A3,B1,A4 form a triangle, with A3 at the centre
c
c      We translate this triangle as a rigid unit

c     ensure we have a valid protein      
      if (iprot.le.0 .or. iprot.gt.NProteins) return

c     get the current and next residue
      ires1 = ires
      ires2 = ires+1

c     ensure that these residues are valid
      if (ires1.gt.NResidues(iprot) .or. ires2.gt.NResidues(iprot)) return
      
c     calculate the translation vector      
      do i=1,3
        ivec(i) = vec(i)*del
      enddo
      
      do i=1,3
c       translate bbatoms 2,3,4 of the first residue
        do j=2,4
          BBCoord(iprot,ires1,j,i) = BBCoord(iprot,ires1,j,i) + ivec(i)
        enddo
c       and then bbatom 1 of the second residue
        BBCoord(iprot,ires2,1,i) = BBCoord(iprot,ires2,1,i) + ivec(i)
      enddo
      
      end
