      subroutine readResTemplateZmatLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c######################################################
c
c     This routine is used to read in zmatrix values 
c     (length of bond/angle/dihedral)
c
c######################################################

      character*(*) line
      integer itmpl,iatm
      integer findTemplateAtom
      double precision readDouble,deg2Rad
      character*4 atmnam
      
c     the format of the line is 'zmat inam bndval angval dihval'
      if (nwords.lt.5) then
        call printLine(WARNING,"Incorrectly formatted residue zmat template line?")
        call printLine(WARNING,line)
        return
      endif
      
c     get the atom name
      atmnam = words(2)(1:4)
      
c     find this atom in the template
      iatm = findTemplateAtom(itmpl,atmnam,0)
      if (iatm.le.0) then
        write(printstring,*) "Cannot find atom ",atmnam," for zmat line"
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
        return
      endif
      
c     save this information with the template atom
      TemplateZmatVal(itmpl,iatm,1) = readDouble(words(3))
      TemplateZmatVal(itmpl,iatm,2) = deg2Rad(readDouble(words(4)))
      TemplateZmatVal(itmpl,iatm,3) = deg2Rad(readDouble(words(5)))
      
      return
      end
