      subroutine readGCsoluteTemplate(gcsnam,FILE,line,changemode,endfile)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c####################################################################
c
c     This controls the reading in of gc-solute templates
c
c     MSB 31/03/09
c
c####################################################################

      integer FILE
      character*(*) line
      character*4 gcsnam
      logical changemode,endfile
      integer itmpl,findGCsoluteTemplate
      
      changemode = .false.
      endfile = .false.
      
c     see if this solvent template has already been loaded
      itmpl = findGCsoluteTemplate(gcsnam)
       if (itmpl.ne.0) then
        write(printstring,*) "Overwriting gc-solute template for solute ",gcsnam
        call printLine(WARNING,printstring)
c       clear the current version of the template
        NGCSTemplateAtms(itmpl) = 0
      else
        itmpl = NGCSTemplates + 1
        if (NGCSTemplates.gt.MAXSVNTEMPLATES) then
          write(printstring,*) "Exceeded maximum number of solvent templates (",
     .                MAXSVNTEMPLATES,")"
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXSVNTEMPLATES and recompile"
          call closeProgram(-1,printstring)
          return
        endif
        NGCSTemplates = itmpl
        endif

      GCSTemplateName(itmpl) = gcsnam
      NGCSTemplateAtms(itmpl) = 0
c     now get on with reading the file
100   continue
        read(FILE,10,err=9700,end=9800) line
10      format(a300)
        
c       skip comments        
        if (line(1:1).eq.'#') goto 100
c       lowercase the line...
        call lowercase(line)
c       split into words
        call split(line)
c       make sure there are at least two words
        if (nwords.lt.2) goto 100
        
        if (words(1).eq.'mode') then
          changemode = .true.
          goto 9900
        else if (words(1).eq.'residue' .or. words(1).eq.'solute'
     .         .or. words(1).eq.'solvent' .or. words(1).eq.'grand') then
         goto 9900
        else if (words(1).eq.'atom') then
          call readGCSTemplateAtomLine(line,itmpl)
        else if (words(1).eq.'info') then
          call readGCSTemplateInfoLine(line,itmpl)
        else
          call printLine(WARNING,"Unrecognised parameter file line (SOLVENT-TEMPLATE)")
          call printLine(WARNING,line)
        endif
c"        
      goto 100
9700  continue
      call printLine(WARNING,"There was an error reading the parameter file (SOLVENT-TEMPLATE)")
9800  continue
      endfile = .true.
9900  continue
c"      
      write(printstring,21) GCSTemplateName(itmpl),NGCSTemplateAtms(itmpl)
      call printLine(DETAIL,printstring)
21    format("GC-solute '",a,"' contains ",i3," atoms.")
      
      return
      end
