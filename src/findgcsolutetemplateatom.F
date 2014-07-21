      integer function findGCsoluteTemplateAtom(inam,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c#########################################################
c
c     This returns the number of atom inam in template
c     itmpl, or 0 if it doesn't exist
c
c     (C) MSB Apr2012 - modded from CW diles
c
c#########################################################

      integer itmpl
      character*4 inam
      integer i,nats
      
      findGCsoluteTemplateAtom = 0
      if (itmpl.le.0 .or. itmpl.gt.NGCSTemplates) return
      
      nats = NGCSTemplateAtms(itmpl)
c      PRINT*, "nats" ,nats
      if (nats.le.0) return
      do i=1,nats
        if (GCSTemplateAtmName(itmpl,i).eq.inam) then
          findGCsoluteTemplateAtom = i
        
          return
        endif
      enddo
      return
      end
