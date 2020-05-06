
      subroutine readResidueOverride(resnm,ipos,FILE,line,changemode,endfile)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'printer.inc'
c###########################################################################
c
c     This routine reads in the override template for protein 
c     residues. Note that the override template can refer to any
c     atoms in the residue, including chain/backbone atoms.
c
c     (C) Christopher Woods, November 2004
c
c###########################################################################

      character*(*) resnm
      integer FILE,ipos
      character*(*) line
      logical changemode,endfile
      integer itmpl,findOverride,iatm
      integer ipar0,ipar1,readInt
      character*6 positionString
      
      changemode = .false.
      endfile = .false.
      
c     ensure that this position is valid
      if (ipos.eq.NOPOS) then
        call printLine(WARNING,"Cannot read a residue override template without a valid chain position!")
        return
      endif
      
c     now see if this residue in this position has already got an override template
      itmpl = findOverride(resnm,ipos)
      
      if (itmpl.eq.0) then
c       we need to add this residue override template
        itmpl = NResOverrides(ipos) + 1
        if (itmpl.gt.MAXTEMPLATES) then
          write(printstring,*) "Exceeded maximum number of override templates (",MAXTEMPLATES,")"
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXTEMPLATES and recompile"
          call closeProgram(-1,printstring)
          return
        endif
        NResOverrides(ipos) = itmpl
c       save the residue name        
        ResOverrideNam(ipos,itmpl) = resnm        
        
c       now create an override template in the override template array
        NOverrides = NOverrides + 1
        if (NOverrides.gt.MAXOVERRIDES) then
          write(printstring,*) "Exceeded maximum number of override templates (",MAXOVERRIDES,")"
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXOVERRIDES and recompile"
          call closeProgram(-1,printstring)
          NOverrides = MAXOVERRIDES
          NResOverrides(ipos) = itmpl-1
          return
        endif
        
        ResOverrideIndex(ipos,itmpl) = NOverrides
        itmpl = NOverrides
      else
        write(printstring,*) "Overwriting the residue override template for ",resnm," in position ",
     .                       positionString(ipos)
        call printLine(WARNING,printstring)
      endif
      
c     initialise this override template
      NOverrideParams(itmpl) = 0
      
c     now read in the atoms...
100   continue
        read(FILE,10,err=9700,end=9800) line
10      format(a300)

        if (line(1:1).eq.'#') goto 100
        
        call lowercase(line)
        call split(line)
        if (nwords.lt.2) goto 100
        
        if (words(1).eq.'mode') then
          changemode = .true.
          goto 9900
        else if (words(1).eq.'residue') then
          goto 9900
        else if (words(1).eq.'atom') then
c         read in an override atom line, format 'atom name ipar0 ipar1' or 'atom name ipar'
          if (nwords.lt.3) then
            call printLine(WARNING,"Cannot interpret residue override atom line!")
            call printLine(WARNING,line)
            goto 100
          else if (nwords.eq.3) then
            ipar0 = readInt(words(3))
            ipar1 = ipar0
          else
            ipar0 = readInt(words(3))
            ipar1 = readInt(words(4))
          endif
          
          iatm = NOverrideParams(itmpl) + 1
          if (iatm.gt.MAXOVERRIDEPARAMS) then
            write(printstring,*) "Exceeded maximum number of override parameters per template (",
     .                            MAXOVERRIDEPARAMS,")"
            call printLine(FATAL,printstring)
            write(printstring,*)"Increase MAXOVERRIDEPARAMS and recompile"
            call closeProgram(-1,printstring)
            goto 100
          endif
          
          NOverrideParams(itmpl) = iatm
          OverrideAtm(itmpl,iatm) = words(2)(1:4)
          OverrideParam(itmpl,iatm,1) = ipar0
          OverrideParam(itmpl,iatm,2) = ipar1
        else
          call printLine(WARNING,"Could not recognise override residue template line")
          call printLine(WARNING,line)
        endif
      goto 100
      
9700  continue
      call printLine(WARNING,"There was an error reading the parameter file (RESIDUE OVERRIDE)")
9800  continue
      endfile = .true.
9900  continue

      write(printstring,20) resnm,positionString(ipos),NOverrideParams(itmpl)
20    format("Residue ",a," in position ",a," has ",i4," override parameters")
      call printLine(DETAIL,printstring)
      
      return
      end
