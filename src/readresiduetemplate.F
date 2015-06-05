      
      subroutine readResidueTemplate(iresnam,FILE,line,changemode,endfile)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'constants.inc'
      include 'enums.inc'
c###############################################################
c
c     This routine is used to read in a residue template
c
c     (C) Christopher Woods  September 2004
c
c###############################################################

      character*4 iresnam
      integer FILE
      character*(*) line
      logical changemode,endfile
      integer itmpl,findResidueTemplate
      integer i,j
      
      changemode = .false.
      endfile = .false.
            
c     see if the residue has been used before
      itmpl = findResidueTemplate(iresnam)
      if (itmpl.ne.0) then
        write(printstring,*) "Overwriting residue template for residue ",iresnam
        call printLine(WARNING,printstring)
c       clear the current version of the template
        NTemplateAtms(itmpl) = 0
        NTemplateBonds(itmpl) = 0
        NTemplateAngles(itmpl) = 0
        NTemplateDihedrals(itmpl) = 0
        NTemplateUreyBradleys(itmpl) = 0
        do i=1,MAXTEMPLATEATMS
          do j=1,3
            TemplateZmatVal(itmpl,i,j) = ZERO
          enddo
        enddo
      else
        NTemplates = NTemplates + 1
        if (NTemplates.gt.MAXTEMPLATES) then
          write(printstring,*) "Exceeded maximum number of residue templates (",
     .                          MAXTEMPLATES,")"
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXTEMPLATES and recompile"
          call closeProgram(-1,printstring)
          NTemplates = MAXTEMPLATES
          return
        endif
        
        itmpl = NTemplates
        TemplateName(itmpl) = iresnam
        
        NTemplateAtms(itmpl) = 0
        NTemplateBonds(itmpl) = 0
        NTemplateAngles(itmpl) = 0
        NTemplateDihedrals(itmpl) = 0
        NTemplateUreyBradleys(itmpl) = 0
        do i=1,MAXTEMPLATEATMS
          do j=1,3
            TemplateZmatVal(itmpl,i,j) = ZERO
          enddo
        enddo
      endif
                      
c     now get to reading the lines
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
        else if (words(1).eq.'residue' .or. words(1).eq.'solute'
     .                   .or. words(1).eq.'solvent' .or. words(1).eq.'chain') then
          goto 9900
        else if (words(1).eq.'atom') then
          call readResTemplateAtomLine(line,itmpl)
        else if (words(1).eq.'bond') then
          call readResTemplateBondLine(line,itmpl)
        else if (words(1).eq.'angle') then
          call readResTemplateAngleLine(line,itmpl)
        else if (words(1).eq.'dihedral') then
          call readResTemplateDihedralLine(line,itmpl)
        else if (words(1).eq.'ureybradley') then
          call readResTemplateUreyBradleyLine(line,itmpl)
        else if (words(1).eq.'zmat') then
          call readResTemplateZmatLine(line,itmpl)
        else if (words(1).eq.'info') then
          call readResTemplateInfoLine(line,itmpl)
        else if (words(1).eq.'backbone') then
          call readResTemplateBackboneLine(line,itmpl)
        else
          call printLine(WARNING,"Unrecognised parameter file line (RESIDUE-TEMPLATE)")
          call printLine(WARNING,line)
        endif

      goto 100
9700  continue
      call printLine(WARNING,"There was an error reading the parameter file (RESIDUE-TEMPLATE)")
9800  continue
      endfile = .true.
9900  continue

      write(printstring,21) iresnam,NTemplateBonds(itmpl),NTemplateAngles(itmpl),
     .   NTemplateDihedrals(itmpl),NTemplateUreyBradleys(itmpl),NTemplateAtms(itmpl)
      call printLine(DETAIL,printstring)
21    format("Residue '",a,"' contains ",i2," bonds, ",i2," angles ",i2," dihedrals ",
     .        i3," Urey-Bradley terms and ",i3," atom types.")
      
      return

      end
