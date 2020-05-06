      subroutine readChainTemplate(inam,FILE,line,changemode,endfile)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'constants.inc'
      include 'enums.inc'
c#####################################################################
c
c     This subroutine reads in the template for a chain
c     backbone.
c
c     (C) Christopher Woods,  October 2004
c
c#####################################################################

      character*(*) inam
      integer FILE
      character*(*) line
      logical changemode,endfile
      integer itmpl,findChainTemplate
      integer i,j,lineLength
      
      changemode = .false.
      endfile = .false.
      
c     see if the chain name has been used before
      itmpl = findChainTemplate(inam)
      
      if (itmpl.ne.0) then
        write(printstring,*) "Overwriting template for backbone chain ",inam
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
c       clear the current version of the template
        NChainAtoms(itmpl) = 0
        NChainBonds(itmpl) = 0
        NChainAngles(itmpl) = 0
        NChainDihedrals(itmpl) = 0
        NChainUreyBradleys(itmpl) = 0
        do i=1,4
          ChainBBAtom(itmpl,i) = "   "
          ChainBBPar(itmpl,i,1) = 0
          ChainBBPar(itmpl,i,2) = 0
        enddo
        do i=1,MAXCHAINATOMS
          do j=1,3
            ChainZmatVal(itmpl,i,j) = ZERO
          enddo
        enddo
      else
c       increment the number of chain templates
        itmpl = NChainTemplates + 1
        if (itmpl.gt.MAXCHAINTEMPLATES) then
          write(printstring,*) "Exceeded maximum number of backbone chain templates (",
     .                          MAXCHAINTEMPLATES,")"
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXCHAINTEMPLATES and recompile"
          call closeProgram(-1,printstring)
          changemode = .true.
          return
        endif
        NChainTemplates = itmpl
        ChainName(itmpl) = inam
c       clear the current version of the template
        NChainAtoms(itmpl) = 0
        NChainBonds(itmpl) = 0
        NChainAngles(itmpl) = 0
        NChainDihedrals(itmpl) = 0
        NChainUreyBradleys(itmpl) = 0
        do i=1,4
          ChainBBAtom(itmpl,i) = "   "
          ChainBBPar(itmpl,i,1) = 0
          ChainBBPar(itmpl,i,2) = 0
        enddo
        do i=1,MAXCHAINATOMS
          do j=1,3
            ChainZmatVal(itmpl,i,j) = ZERO
          enddo
        enddo
      endif
      
c     read in the template
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
        else if (words(1).eq.'bbatom') then
          call readChainTemplateBBAtomLine(line,itmpl)
        else if (words(1).eq.'atom') then
          call readChainTemplateAtomLine(line,itmpl)
        else if (words(1).eq.'bond') then
          call readChainTemplateBondLine(line,itmpl)
        else if (words(1).eq.'angle') then
          call readChainTemplateAngleLine(line,itmpl)
        else if (words(1).eq.'dihedral') then
          call readChainTemplateDihedralLine(line,itmpl)
        else if (words(1).eq.'ureybradley') then
          call readChainTemplateUBLine(line,itmpl)
        else if (words(1).eq.'zmat') then
          call readChainTemplateZmatLine(line,itmpl)
        else
          call printLine(WARNING,"Unrecognised parameter file line (CHAIN-TEMPLATE)")
          call printLine(WARNING,line)
        endif

      goto 100
9700  continue
      call printLine(WARNING,"There was an error reading the parameter file (CHAIN-TEMPLATE)")
9800  continue
      endfile = .true.
      
9900  continue
      
      i = lineLength(inam)
      write(printstring,21) inam(1:i),ChainBBAtom(itmpl,1),ChainBBAtom(itmpl,2),
     .   ChainBBAtom(itmpl,3),ChainBBAtom(itmpl,4),NChainBonds(itmpl),NChainAngles(itmpl),
     .   NChainDihedrals(itmpl),NChainUreyBradleys(itmpl),NChainAtoms(itmpl)
      call printLine(DETAIL,printstring)
21    format("Chain '",a,"' has backbone atoms (",a,"-",a,"-",a,"-",a,"), contains ",
     .        i2," bonds, ",i2," angles ",i2," dihedrals ",
     .        i3," Urey-Bradley terms and ",i3," atom types.")
      
      return

      end
