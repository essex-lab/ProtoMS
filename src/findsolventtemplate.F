      integer function findSolventTemplate(inam)
      include 'dimensions.inc'
      include 'templates.inc'
c#############################################################
c
c     This function returns the template number of solvent inam
c     (or 0 if it doesn't exist)
c
c     (C) Christopher Woods   5 March 2002
c
c#############################################################

      character*4 inam
      integer i
      
      findSolventTemplate = 0
      
      if (NSvnTemplates.le.0) return
      
      do i=1,NSvnTemplates
        if (SvnTemplateName(i).eq.inam) then
          findSolventTemplate = i
          return
        endif
      enddo

      return
      end
