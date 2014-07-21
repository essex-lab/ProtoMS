      integer function findTemplate(inam)
      include 'dimensions.inc'
      include 'templates.inc'
c#######################################################
c
c     This function will attempt to find the template
c     that belongs to residue with resname
c
c     (C) Christopher Woods, October 2004
c
c#######################################################

      character*4 inam
      integer i

c     loop over all templates, and see if we have
c     a match...
      findTemplate = 0
      if (NTemplates.eq.0) return

      do i=1,NTemplates
        if (TemplateName(i).eq.inam) then
          findTemplate = i
          return
        endif
      enddo

      return
      end
