      integer function findGCsoluteTemplate(inam)
      include 'dimensions.inc'
      include 'templates.inc'
c#############################################################
c
c     This function returns the template number of gc-solute inam
c     (or 0 if it doesn't exist)
c
c     (C) MSB 1 April 2009
c
c#############################################################

      character*4 inam
      integer i
      
      findGCsoluteTemplate = 0
      
      if (NGCSTemplates.le.0) return
      
      do i=1,NGCSTemplates
        if (GCSTemplateName(i).eq.inam) then
          findGCsoluteTemplate = i
          return
        endif
      enddo

      return
      end
