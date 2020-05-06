      integer function getChainTemplate(itmpl,ipos)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'templates.inc'
c###########################################################
c
c     This function returns the chain template for 
c     residue template itmpl in chain position ipos
c
c     0 is returned if there is no chain template 
c     for this residue in this position
c
c     (C) Christopher Woods
c
c##########################################################

      integer itmpl
      integer ipos,findChainTemplate
      
      getChainTemplate = 0
      
      if (ipos.eq.SINGLE .or. ipos.eq.FIRST .or.
     .    ipos.eq.MIDDLE .or. ipos.eq.LAST) then
        if (TemplateResChain(itmpl,ipos).eq.'none') return
        
        getChainTemplate = findChainTemplate(TemplateResChain(itmpl,ipos))
        return
      else
        return
      endif
      
      end
