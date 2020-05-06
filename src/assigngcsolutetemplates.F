      subroutine assigngcsoluteTemplates()
      include 'dimensions.inc'
      include 'solvents.inc'
c#####################################################
c
c     This routine assigns all of the gc-solute molecule
c     templates!
c
c     MSB 1st April 09
c
c#####################################################

      integer i

      if (NGCsolutes.le.0) return
      
      do i=1,NGCsolutes
        call assigngcsoluteTemplate(i)
      
       enddo
      
      return
      end
