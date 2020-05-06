      subroutine readSoluteTemplate(solnm,FILE,line,changemode,endfile)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
      include 'constants.inc'
c###############################################################
c
c     This is a master routine that coordinates the reading in
c     of a solute template
c
c     (C) Christopher Woods,  September 2004
c
c     Modified by Samuel Genheden 2012-2013
c
c###############################################################

      character*(*) solnm
      integer FILE
      character*(*) line
      logical changemode,endfile
      integer itmpl,i
      integer lineLength,lgth
      integer findMoleculeTemplate
      double precision readDouble
      
      changemode = .false.
      endfile = .false.
      
c     see if the solute name has been used before
      itmpl = findMoleculeTemplate(solnm)
      if (itmpl.ne.0) then
        write(printstring,*) "Overwriting solute template for solute ",solnm
        call printLine(WARNING,printstring)
c       clear the current version of the solute template
        NMolTemplateBonds(itmpl) = 0
        NMolTemplateAngles(itmpl) = 0
        NMolTemplateDihedrals(itmpl) = 0
        NMolTemplateUreyBradleys(itmpl) = 0
        NMolTemplateAtms(itmpl) = 0      
      else
        itmpl = NMolTemplates + 1
        if (itmpl.gt.MAXMOLTEMPLATES) then
          write(printstring,*) "Exceeded maximum number of molecule templates (",
     .                         MAXMOLTEMPLATES,")"
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXMOLTEMPLATES and recompile"
          call closeProgram(-1,printstring)
          return
        endif
        NMolTemplates = itmpl
      endif
      
      MolTemplateName(itmpl) = solnm
      
c     Zeroing JAWS variable in case it is not in the template file
      MolTemplateTDel(itmpl) = ZERO
      MolTemplateJCorr(itmpl) = ZERO
      MolTemplateJPmfOrd(itmpl) = 0

c     now get on with reading the lines...
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
     .               .or. words(1).eq.'solvent') then
          goto 9900
        else if (words(1).eq.'atom') then
          call readSolTemplateAtomLine(line,itmpl)
        else if (words(1).eq.'bond') then
          call readSolTemplateBondLine(line,itmpl)
        else if (words(1).eq.'angle') then
          call readSolTemplateAngleLine(line,itmpl)
        else if (words(1).eq.'dihedral') then
          call readSolTemplateDihedralLine(line,itmpl)
        else if (words(1).eq.'ureybradley') then
          call readSolTemplateUreyBradleyLine(line,itmpl)
        else if (words(1).eq.'variable') then
          call readSolTemplateVariableLine(line,itmpl)
        else if (words(1).eq.'coord') then
          call readSolTemplateCoordLine(line,itmpl)
        else if (words(1).eq.'info') then
          call readSolTemplateInfoLine(line,itmpl)
!       SG, Read JAWS-specific variables
        else if (words(1).eq.'jtheta') then
          MolTemplateTDel(itmpl) = readDouble(words(2),0.0)
        else if (words(1).eq.'jpmf') then
          MolTemplateJPmfOrd(itmpl) = nwords-1
          if (MolTemplateJPmfOrd(itmpl).gt.MAXJPMFORDER) MolTemplateJPmfOrd(itmpl) = MAXJPMFORDER
          do i=2,nwords
            MolTemplateJPmfCoeff(itmpl,i-1) = readDouble(words(i),0.0)
          enddo
        else if (words(1).eq.'jcorr') then
          MolTemplateJCorr(itmpl) = readDouble(words(2),0.0)
        else
          call printLine(WARNING,"Unrecognised parameter file line (SOLUTE-TEMPLATE)")
          call printLine(WARNING,line)
        endif
        
      goto 100
9700  continue
      call printLine(WARNING,"There was an error reading the parameter file (SOLUTE-TEMPLATE)")
9800  continue
      endfile = .true.
9900  continue

      lgth = lineLength(solnm)

      write(printstring,21) solnm(1:lgth),NMolTemplateBonds(itmpl),NMolTemplateAngles(itmpl),
     .  NMolTemplateDihedrals(itmpl),NMolTemplateUreyBradleys(itmpl),
     .  NMolTemplateAtms(itmpl)
21    format("Solute '",a,"', contains ",i5," bonds, ",i5," angles ",i5," dihedrals ",i5,
     .          " Urey-Bradley terms and ",i5," atom types.")
      call printLine(DETAIL,printstring)
      return
          
      end
