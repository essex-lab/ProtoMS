      integer function findTemplateBond(at1,at2,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c#####################################################
c
C     Same as findBond, but on a4 types, and on template
c     array...
c
c     (C) Christopher Woods  23 January 2002
c
c######################################################

      character*4 at1,at2
      integer i,itmpl
      
      findTemplateBond = 0
      if (NTemplateBonds(itmpl).le.0) return

      call sortTemplateBond(at1,at2)
      
      do i=1,NTemplateBonds(itmpl)
        if (at1.eq.TemplateBndName(itmpl,i,1)) then
          if (at2.eq.TemplateBndName(itmpl,i,2)) then
            findTemplateBond = i
            return
          endif
        endif
      enddo

      return
      end
