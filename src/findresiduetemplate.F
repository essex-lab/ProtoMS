      integer function findResidueTemplate(inam)
      include 'dimensions.inc'
      include 'templates.inc'
c###############################################
c
c     This function returns the index of the 
c     template residue inam. It returns 0 if
c     it doesn't exist in the template list
c
c     (C) Christopher Woods   23 January 2002
c
c###############################################

      character*4 inam
      integer i

      findResidueTemplate = 0
      
      if (NTemplates.le.0) return
      
      do i=1,NTemplates
        if (TemplateName(i).eq.inam) then
           findResidueTemplate = i
           return
        endif
      enddo

      return
      
      end
