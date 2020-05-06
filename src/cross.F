      subroutine cross(vec1,vec2,norm)
      implicit none
c##################################################################
c
c     This routine calculates the normal (norm) for vectors
c     vec1 and vec2 from the cross product - the normal is
c     normalised
c
c     (C) Christopher Woods         21 August 2002
c
c##################################################################

      double precision norm(3),vec1(3),vec2(3)
      
      norm(1) = (vec1(2)*vec2(3)) - (vec1(3)*vec2(2))
      norm(2) = (vec1(3)*vec2(1)) - (vec1(1)*vec2(3))
      norm(3) = (vec1(1)*vec2(2)) - (vec1(2)*vec2(1))
      
      call vecnormalise(norm)
      
      return
      end
