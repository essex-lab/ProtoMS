      integer function findTemplateAngle(at1,at2,at3,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c#####################################################
c
c     Same as findangle, but for templates
c
c     (C) Christopher Woods  23 January 2002
c
c######################################################

      character*4 at1,at2,at3
      integer i,itmpl

      findTemplateAngle = 0
      if (NTemplateAngles(itmpl).le.0) return

      call sortTemplateAngle(at1,at2,at3)
      
      do i=1,NTemplateAngles(itmpl)
        if (at1.eq.TemplateAngName(itmpl,i,1)) then
          if (at2.eq.TemplateAngName(itmpl,i,2)) then
            if (at3.eq.TemplateAngName(itmpl,i,3)) then
              findTemplateAngle = i
              return
            endif
          endif
        endif
      enddo

      return
      end
