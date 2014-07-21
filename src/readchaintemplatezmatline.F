      subroutine readChainTemplateZmatLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'splitter.inc'
c######################################################
c
c     This routine is used to read in zmatrix values 
c     (length of bond/angle/dihedral)
c
c######################################################

      character*(*) line
      integer itmpl,iatm
      integer findChainAtom
      double precision readDouble,deg2Rad
      character*4 atmnam
      
c     the format of the line is 'zmat inam bndval angval dihval'
      if (nwords.lt.5) then
        call printLine(WARNING,"Incorrectly formatted chain zmat template line?")
        call printLine(WARNING,line)
        return
      endif
      
c     get the atom name
      atmnam = words(2)(1:4)
      
c     find this atom in the template
      iatm = findChainAtom(itmpl,atmnam)
      if (iatm.eq.0) then
        write(printstring,*) "Cannot find atom ",atmnam," for zmat line"
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
        return
      else if (iatm.lt.0) then
        write(printstring,*) "Cannot have a zmat line for bbatom ",atmnam
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
        return
      endif
      
c     save this information with the template atom - convert angles to radians
      ChainZmatVal(itmpl,iatm,1) = readDouble(words(3))      
      ChainZmatVal(itmpl,iatm,2) = deg2Rad(readDouble(words(4)))
      ChainZmatVal(itmpl,iatm,3) = deg2Rad(readDouble(words(5)))
      
      return
      end
