      
      subroutine readSvnTemplateAtomLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c###############################################################
c
c     Read in a solvent atom template line
c
c     (C) Christopher Woods  September 2004
c
c###############################################################

      character*(*) line
      integer itmpl
      character*4 inam
      integer ipar0,ipar1
      integer iatm
      integer readInt,findSolventTemplateAtom
      
c     the format is "atom name ipar0 ipar1"
      if (nwords.ge.4) then
        ipar0 = readInt(words(3))
        ipar1 = readInt(words(4))
      else if (nwords.eq.3) then
c       format is "atom name ipar"  ipar0=ipar1=ipar
        ipar0 = readInt(words(3))
        ipar1 = ipar0
      else if (nwords.eq.2) then
c       format is "atom name" ipar0=ipar1=0
        ipar0 = 0
        ipar1 = 0
      else
        call printLine(WARNING,"Incorrectly formatted solvent atom template line?")
        call printLine(WARNING,line)
        return
      endif
      
c     read the data itself
      inam = words(2)(1:4)

c     make sure that clj parameters below 0 are interpreted as dummies
      if (ipar0.lt.0 .or. ipar0.gt.MAXCLJ) ipar0 = 0
      if (ipar1.lt.0 .or. ipar1.gt.MAXCLJ) ipar1 = 0
              
c     make sure that this atom doesn't already exist
      iatm = findSolventTemplateAtom(inam,itmpl)
      if (iatm.ne.0) then
        write(printstring,*) "Overwriting solvent template atom data for ",iatm
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
      else
c       increment the number of atoms
        iatm = NSvnTemplateAtms(itmpl) + 1
        if (iatm.gt.MAXSVNTEMPLATEATMS) then
          write(printstring,*) "Exceeded maximum number of solvent template atoms (",
     .                  MAXSVNTEMPLATEATMS,")"
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXSVNTEMPLATEATMS and recompile"
          call closeProgram(-1,printstring)
          return
        endif
        NSvnTemplateAtms(itmpl) = iatm
      endif
      
      SvnTemplateAtmName(itmpl,iatm) = inam
      SvnTemplateAtmCLJ(itmpl,iatm,1) = ipar0
      SvnTemplateAtmCLJ(itmpl,iatm,2) = ipar1
      
      return
      end
