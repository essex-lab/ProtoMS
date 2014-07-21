      subroutine assignSolventTemplates()
      include 'dimensions.inc'
      include 'solvents.inc'
c#####################################################
c
c     This routine assigns all of the solvent molecule
c     templates!
c
c     (C) Christopher Woods  5 March 2002
c
c#####################################################

      integer i

      if (NSolvents.le.0) return
      
      do i=1,NSolvents
        call assignSolventTemplate(i)
      enddo
      
      return
      end
