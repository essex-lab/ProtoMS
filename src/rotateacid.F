      subroutine rotateAcid(iprot,ires,ori,vec,del)
      include 'dimensions.inc'
      include 'proteins.inc'
c####################################################
c
c     This function translates the backbone
c     atoms of residue ires in protein iprot
c     by del radians around the vector vec 
c     with the origin at ori
c
c     (C) Christopher Woods, October 2004
c
c####################################################
      
      double precision ori(3),vec(3)
      double precision ivec(3)
      double precision del,ilgth
      double precision vecLength
      integer iprot,ires,ires1,i,j

c     backbone atoms form a triangle with next residue...
c     (residue A, bb atoms 1-4, residue B, bb atoms 1-4
c      
c              A4          B4
c              |           |
c      A1--A2--A3--B1--B2--B3
c
c     So A2,A3,B1,A4 form a triangle, with A3 at the centre, e.g
c
c             A4
c             |
c             A3
c            /  \
c          A2 .. B1
c
c     We rotate this triangle as a rigid unit around the point 'ori' 
c     (here marked as a '..') around the vector ivec (normally
c     parallel to vector A2-B1, the parallel to A3-A4, or
c     perpendicular to the plane of the triangle)
      
      if (iprot.le.0 .or. iprot.gt.NProteins) return
      if (ires.le.0 .or. ires.ge.NResidues(iprot)) return

c     do each atom...
      do j=1,4
        if (j.eq.1) then
          ires1 = ires+1
        else
          ires1 = ires
        endif
      
c       get the vector from the origin to the atom
        do i=1,3
          ivec(i) = BBCoord(iprot,ires1,j,i) - ori(i)
        enddo

c       save the length of this vector      
        ilgth = vecLength(ivec)

c       rotate this vector
        call rotateVector(ivec,vec,del)

c       reset the length of the vectors - this prevents
c       precision errors from propagating
        call setLength(ilgth,ivec)
      
c       set the new coordinates of the atom
        do i=1,3
          BBCoord(iprot,ires1,j,i) = ori(i) + ivec(i)
        enddo
      enddo
              
      return
      end
