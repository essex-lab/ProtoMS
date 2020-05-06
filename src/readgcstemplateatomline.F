      
      subroutine readGCSTemplateAtomLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c###############################################################
c
c     Read in a gc-solute atom template line
c
c     MSB 31/03/09
c
c###############################################################

      character*(*) line
      integer itmpl
      character*4 inam
      integer ipar0,ipar1
      integer iatm
      integer readInt,findGCsoluteTemplateAtom
      
c     the format is "atom name ipar0 ipar1"
      if (nwords.ge.4) then
        ipar0 = readInt(words(3))
c        PRINT*, "ipar0" ,ipar0
        ipar1 = readInt(words(4))
      else if (nwords.eq.3) then
c      format is "atom name ipar"  ipar0=ipar1=ipar
        ipar0 = readInt(words(3))
        ipar1 = ipar0
      else if (nwords.eq.2) then
c         PRINT*, "TWO"
c       format is "atom name" ipar0=ipar1=0
        ipar0 = 0
        ipar1 = 0
      else
        call printLine(WARNING,"Incorrectly formatted GCsolute atom template line?")
        call printLine(WARNING,line)
        return
      endif
      
c     read the data itself
      inam = words(2)(1:4)
     
     
c     make sure that clj parameters below 0 are interpreted as dummies
      if (ipar0.lt.0 .or. ipar0.gt.MAXCLJ) ipar0 = 0
      if (ipar1.lt.0 .or. ipar1.gt.MAXCLJ) ipar1 = 0
c"
       
              
c     make sure that this atom doesn't already exist
      iatm = findGCsoluteTemplateAtom(inam,itmpl)
       if (iatm.ne.0) then
        write(printstring,*) "Overwriting GCsolute template atom data for ",iatm
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
      
       else
c       increment the number of atoms
       
        iatm = NGCSTemplateAtms(itmpl) + 1
        
        if (iatm.gt.100) then
          write(printstring,*) "Exceeded maximum number of GCsolute template atoms (",
     .                  100,")"
          call printLine(WARNING,printstring)
          call printLine(WARNING,"Increase MAXSOLTEMPLATEATMS and recompile")
          return
        endif
        NGCSTemplateAtms(itmpl) = iatm
c        PRINT*,"OUT" ,NGCSTemplateAtms(itmpl)
      endif
      
      GCSTemplateAtmName(itmpl,iatm) = inam
      GCSTemplateAtmCLJ(itmpl,iatm,1) = ipar0
      GCSTemplateAtmCLJ(itmpl,iatm,2) = ipar1
c      PRINT*, "debug" ,GCSTemplateAtmName(itmpl,iatm) , GCSTemplateAtmCLJ(itmpl,iatm,1), GCSTemplateAtmCLJ(itmpl,iatm,2)
c      PRINT*, "test" ,GCSTemplateAtmCLJ(itmpl,iatm,1)
c      PRINT*, "test" ,GCSTemplateAtmName(itmpl,iatm)

      return
      end
