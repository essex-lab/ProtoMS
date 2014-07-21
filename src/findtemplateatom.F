      integer function findTemplateAtom(itmpl,inam,ipos)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'enums.inc'
c#######################################################
c
c     This function finds and returns the index of
c     the atom called 'inam' in the residue template
c     'itmpl'. 0 is returned if the atom does not
c     exist. 
c
c     This function amalgamates the residue template
c     (itmpl) with the chain template for chain
c     position ipos (FIRST,MIDDLE,LAST,SINGLE). If this
c     atom is a bbatom, then -1 to -4 will be 
c     returned (depending on the BBAtom type)
c
c     (C) Christopher Woods, October 2004
c
c######################################################

      integer itmpl
      character*4 inam
      integer ipos
      integer ichn
      integer natms,i
      integer getChainTemplate,findChainAtom

      findTemplateAtom = 0
      
      if (itmpl.eq.0) then
        call printLine(WARNING,"Cannot search template 0 - this is certainly a bug!")
        return
      endif
      
c     search the template for this atom
      natms = NTemplateAtms(itmpl)
      
      do i=1,natms
        if (TemplateAtmName(itmpl,i,1).eq.inam) then
          findTemplateAtom = i
          return
        endif
      enddo
      
c     now search the backbone atoms - get the chain template for this
c     residue in this position
      if (ipos.ne.0) then
        ichn = getChainTemplate(itmpl,ipos)
        if (ichn.eq.0) then
          write(printstring,*) "Cannot find chain template for residue ",inam,
     .                         "in position ",ipos,". THIS IS A BUG!"
          call printLine(WARNING,printstring)
        else
              
          i = findChainAtom(ichn,inam)
          if (i.gt.0) then
c           add this number onto the number of atoms in the main template
            findTemplateAtom = i + natms
          else
            findTemplateAtom = i
          endif
        endif
      endif
        
      return
      end
