      
      integer function findMoleculeTemplate(inam)
      include 'dimensions.inc'
      include 'templates.inc'
      include 'printer.inc'
c########################################################
c
c     This function returns the template number of molecule
c     inam (or 0 if it doesn't exist)
c
c     (C) Christopher Woods   18 February 2002
c
c########################################################

      character*(*) inam
      integer lineLength,lgth
      integer i

      findMoleculeTemplate = 0
      
      if (NMolTemplates.le.0) return
      
      lgth = lineLength(inam)
      
      do i=1,NMolTemplates
        if (lineLength(MolTemplateName(i)).eq.lgth) then
          if (MolTemplateName(i)(1:lgth).eq.inam(1:lgth)) then
            findMoleculeTemplate = i
            return
          endif
        endif
      enddo

      write(printstring,*) "Failed to find molecular template when searching for '",
     .                     inam(1:lgth),"'"
      call printLine(DETAIL,printstring)
      
      
      return
      end
