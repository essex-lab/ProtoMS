      integer function findTemplateUreyBradley(at1,at2,at3,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c#####################################################
c
C     Same as findUreyBradley, but on a4 types, and on template
c     array...
c
c     (C) Christopher Woods
c
c######################################################

      character*4 at1,at2,at3
      integer i,itmpl
      
      findTemplateUreyBradley = 0
      if (NTemplateUreyBradleys(itmpl).le.0) return

      call sortTemplateUreyBradley(at1,at2,at3)
      
      do 100 i=1,NTemplateUreyBradleys(itmpl)
        if (at1.eq.TemplateUByName(itmpl,i,1)) then
          if (at2.eq.TemplateUByName(itmpl,i,2)) then
            if (at3.eq.TemplateUByName(itmpl,i,3)) then
              findTemplateUreyBradley = i
              return
            endif
          endif
        endif
100   continue

      return
      end
