      subroutine readSolventTemplate(svnnam,FILE,line,changemode,endfile)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c####################################################################
c
c     This controls the reading in of solvent templates
c
c     (C) Christopher Woods   September 2004
c
c####################################################################

      integer FILE
      character*(*) line
      character*4 svnnam
      logical changemode,endfile
      integer itmpl,findSolventTemplate
      
      changemode = .false.
      endfile = .false.
      
c     see if this solvent template has already been loaded
      itmpl = findSolventTemplate(svnnam)
      
      if (itmpl.ne.0) then
        write(printstring,*) "Overwriting solvent template for solvent ",svnnam
        call printLine(WARNING,printstring)
c       clear the current version of the template
        NSvnTemplateAtms(itmpl) = 0
      else
        itmpl = NSvnTemplates + 1
        if (NSvnTemplates.gt.MAXSVNTEMPLATES) then
          write(printstring,*) "Exceeded maximum number of solvent templates (",
     .                MAXSVNTEMPLATES,")"
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXSVNTEMPLATES and recompile"
          call closeProgram(-1,printstring)
          return
        endif
        NSvnTemplates = itmpl
      endif

      SvnTemplateName(itmpl) = svnnam
c     Default is no optimization
      SvnTemplateOptimization(itmpl) = OPTNONE

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
     .                 .or. words(1).eq.'solvent') then
          goto 9900
        else if (words(1).eq.'atom') then
          call readSvnTemplateAtomLine(line,itmpl)
        else if (words(1).eq.'info') then
          call readSvnTemplateInfoLine(line,itmpl)
        else if (words(1).eq.'opt') then
          if (words(2)(1:5).eq.'tip4p') then
            call printLine(INFO,"Found solvent template with TIP4p optimization")
            SvnTemplateOptimization(itmpl) = OPTTIP4P
          else if (words(2)(1:5).eq.'tip3p') then
            call printLine(INFO,"Found solvent template with TIP3p optimization")            
            SvnTemplateOptimization(itmpl) = OPTTIP3P
          else
            SvnTemplateOptimization(itmpl) = OPTNONE
            call printLine(INFO,"Found solvent template without optimization")          
          endif
        else if (words(1).eq.'coord') then
          call readSvnTemplateCoordLine(line,itmpl)
        else
          call printLine(WARNING,"Unrecognised parameter file line (SOLVENT-TEMPLATE)")
          call printLine(WARNING,line)
        endif
        
      goto 100
9700  continue
      call printLine(WARNING,"There was an error reading the parameter file (SOLVENT-TEMPLATE)")
9800  continue
      endfile = .true.
9900  continue
      
      write(printstring,21) SvnTemplateName(itmpl),NSvnTemplateAtms(itmpl)
      call printLine(DETAIL,printstring)
21    format("Solvent '",a,"' contains ",i3," atoms.")
      
      return
      end
