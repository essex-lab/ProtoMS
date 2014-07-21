      subroutine readResTemplateBackboneLine(line,itmpl)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'templates.inc'
      include 'splitter.inc'
c############################################################
c
c     This routine reads in which backbones are associated with
c     this residue template
c
c     (C) Christopher Woods
c
c############################################################

      character*(*) line
      integer itmpl
      integer i,itmp
      
      TemplateResChain(itmpl,FIRST) = "none"
      TemplateResChain(itmpl,MIDDLE) = "none"
      TemplateResChain(itmpl,LAST) = "none"
      TemplateResChain(itmpl,SINGLE) = "none"

c     format is "backbone (first name) (middle name) (last name) (single name)"
      do i=2,nwords-1
        if (words(i).eq.'first') then
          itmp = FIRST
        else if (words(i).eq.'middle') then
          itmp = MIDDLE
        else if (words(i).eq.'last') then
          itmp = LAST
        else if (words(i).eq.'single') then
          itmp = SINGLE
        else
          itmp = 0
        endif
        
        if (itmp.ne.0) TemplateResChain(itmpl,itmp) = words(i+1)(1:20)
      enddo
      return
      end
