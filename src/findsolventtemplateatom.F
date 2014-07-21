      integer function findSolventTemplateAtom(inam,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c#########################################################
c
c     This returns the number of atom inam in template
c     itmpl, or 0 if it doesn't exist
c
c     (C) Christopher Woods  5 March 2002
c
c#########################################################

      integer itmpl
      character*4 inam
      integer i,nats
      
      findSolventTemplateAtom = 0
      if (itmpl.le.0 .or. itmpl.gt.NSvnTemplates) return
      
      nats = NSvnTemplateAtms(itmpl)
      if (nats.le.0) return
      do i=1,nats
        if (SvnTemplateAtmName(itmpl,i).eq.inam) then
          findSolventTemplateAtom = i
          return
        endif
      enddo
      return
      end
