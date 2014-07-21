      integer function findChainTemplate(inam)
      include 'dimensions.inc'
      include 'templates.inc'
c#######################################################
c
c     This function will attempt to find the template
c     that belongs to backbone chain with name 'inam'
c
c     (C) Christopher Woods  October 2004
c
c#######################################################

      character*(*) inam
      integer lineLength,ilgth,i

c     loop over all templates, and see if we have
c     a match...
      findChainTemplate = 0
      if (NChainTemplates.eq.0) return

      ilgth = lineLength(inam)
      
      do i=1,NChainTemplates
        if (ChainName(i)(1:ilgth).eq.inam(1:ilgth)) then
          findChainTemplate = i
          return
        endif
      enddo

      return
      end
