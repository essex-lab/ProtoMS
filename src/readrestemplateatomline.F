      
      subroutine readResTemplateAtomLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c###################################################
c
c     This routine interprets a set of words that 
c     contains a residue atom template
c
c     (C) Christopher Woods
c
c###################################################

      character*(*) line
      integer itmpl
      integer findTemplateAtom,readInt
      integer nat,ipar0,ipar1,istrt
      integer ipar2,ipar3,ipar4,ipar5
      character*4 at1
            
c     we only need to get the atom name and clj parameter
c     as the rest of the information in the template is used to 
c     build the prozmat file (and is contained in the resulting
c     prozmat file)

c       JM Mod January 2005
c       AMBER can have different CLJ parameters for sidechain atoms depending
c       on the position of the backbone, so we need 6 par numbers, the order 
c       is ''middle'', ''first'', ''last''. The strategy adopted is 
c       to extend the size
c       of TemplateAtmCLJ and fill it even if some values won't be used
c       The correct selection of the appropriate parameter is made in 
c       assignresiduetemplate because we need to know the position of the 
c       residue.

      if (nwords .eq. 11) then
c       AMBER
c       format is "atom name ipar0 ipar1 ipar2 ipar3 ipar4 ipar5 bndatm angatm dihatm"
        ipar0 = readInt(words(3))
        ipar1 = readInt(words(4))
        ipar2 = readInt(words(5))
        ipar3 = readInt(words(6))
        ipar4 = readInt(words(7))
        ipar5 = readInt(words(8))
        istrt = 9
      else if (nwords.eq.7) then
c       format is "atom name ipar0 ipar1 bndatm angatm dihatm"
        ipar0 = readInt(words(3))
        ipar1 = readInt(words(4))
        ipar2 = ipar0
        ipar3 = ipar1
        ipar4 = ipar0
        ipar5 = ipar1
        istrt = 5
      else if (nwords.eq.6) then
c       format is "atom name ipar bndatm angatm dihatm"
        ipar0 = readInt(words(3))
        ipar1 = ipar0
        ipar2 = ipar0
        ipar3 = ipar1
        ipar4 = ipar0
        ipar5 = ipar1
        istrt = 4
      else if (nwords.eq.5) then
c       format is "atom name bndatm angatm dihatm", ipar0 and ipar1 implicitly = 0
        ipar0 = 0
        ipar1 = 0
        ipar2 = ipar0
        ipar3 = ipar1
        ipar4 = ipar0
        ipar5 = ipar1
        istrt = 3
      else
        call printLine(WARNING,"Incorrectly formatted residue atom template line?")
        call printLine(WARNING,line)
        return
      endif

      at1 = words(2)(1:4)
      
      nat = findTemplateAtom(itmpl,at1,0)
      if (nat.eq.0) then
        nat = NTemplateAtms(itmpl) + 1
        if (nat .gt. MAXTEMPLATEATMS) then
          write(printstring,*) "Exceeded maximum number of template atoms (",MAXTEMPLATEATMS,")"
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXTEMPLATEATMS and recompile"
          call closeProgram(-1,printstring)
          return
        endif
        NTemplateAtms(itmpl) = nat      
      else
        write(printstring,*) "Overwriting residue template atom ",at1
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
      endif
      
c     ensure that values less than 0 are interpreted as dummies
      if (ipar0.lt.0 .or. ipar0.gt.MAXCLJ) ipar0 = 0
      if (ipar1.lt.0 .or. ipar1.gt.MAXCLJ) ipar1 = 0
      if (ipar2.lt.0 .or. ipar2.gt.MAXCLJ) ipar2 = 0
      if (ipar3.lt.0 .or. ipar3.gt.MAXCLJ) ipar3 = 0
      if (ipar4.lt.0 .or. ipar4.gt.MAXCLJ) ipar4 = 0
      if (ipar5.lt.0 .or. ipar5.gt.MAXCLJ) ipar5 = 0
      
      TemplateAtmName(itmpl,nat,1) = at1
      TemplateAtmCLJ(itmpl,nat,1) = ipar0
      TemplateAtmCLJ(itmpl,nat,2) = ipar1
      TemplateAtmCLJ(itmpl,nat,3) = ipar2
      TemplateAtmCLJ(itmpl,nat,4) = ipar3
      TemplateAtmCLJ(itmpl,nat,5) = ipar4
      TemplateAtmCLJ(itmpl,nat,6) = ipar5
      TemplateAtmName(itmpl,nat,2) = words(istrt)(1:4)
      TemplateAtmName(itmpl,nat,3) = words(istrt+1)(1:4)
      TemplateAtmName(itmpl,nat,4) = words(istrt+2)(1:4)
      
      return
      end
