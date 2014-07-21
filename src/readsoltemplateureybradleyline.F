      subroutine readSolTemplateUreyBradleyLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c###########################################################
c
c     This routine reads in the Urey-Bradley terms
c     for a molecule template
c
c     (C) Christopher Woods
c
c###########################################################

      character*(*) line
      integer itmpl
      character*4 iatm1,ires1,iatm2,ires2,iatm3,ires3
      logical flex,dummy
      double precision delta
      integer ipar0,ipar1
      integer findMolTemplateUreyBradley
      integer iuby

c     format is;
c     'ureybradley iatm1 ires1 iatm2 ires2 iatm3 ires3 (dummy) (param ipar0 ipar1)'

      if (nwords.lt.7) then
        call printLine(WARNING,"Incorrectly formatted solute template Urey-Bradley line?")
        call printLine(WARNING,line)
        return
      endif
      
c     read in the data
      iatm1 = words(2)(1:4)
      ires1 = words(3)(1:4)
      iatm2 = words(4)(1:4)
      ires2 = words(5)(1:4)
      iatm3 = words(6)(1:4)
      ires3 = words(7)(1:4)
      
      call getDOFInfo(dummy,flex,delta,ipar0,ipar1)

      if (flex) then
        call printLine(WARNING,"Cannot have flexible Urey-Bradley terms!")
        call printLine(WARNING,line)
      endif
      
c     see if this uby already exists
      iuby = findMolTemplateUreyBradley(iatm1,ires1,iatm2,ires2,iatm3,ires3,itmpl)
      if (iuby.ne.0) then
        write(printstring,*) "Overwriting information for solute Urey-Bradley (",
     .      iatm1,",",ires1,")-(",iatm2,",",ires2,")-(",iatm3,",",ires3,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
      else
        iuby = NMolTemplateUreyBradleys(itmpl) + 1
        if (iuby.gt.MAXMOLTEMPLATEUREYBRADLEYS) then
          write(printstring,*) "Exceeded maximum number of molecule template Urey-Bradleys (",
     .                            MAXTEMPLATEUREYBRADLEYS,")"
          call printLine(WARNING,printstring)
          call printLine(WARNING,"Increase MAXTEMPLATEUREYBRADLEYS and recompile")
          return
        endif
        NMolTemplateUreyBradleys(itmpl) = iuby
      endif
       
      MolTemplateUByName(itmpl,iuby,1) = iatm1
      MolTemplateUByRes(itmpl,iuby,1) = ires1
      MolTemplateUByName(itmpl,iuby,2) = iatm2
      MolTemplateUByRes(itmpl,iuby,2) = ires2
      MolTemplateUByName(itmpl,iuby,3) = iatm3
      MolTemplateUByRes(itmpl,iuby,3) = ires3
      MolTemplateUByPar(itmpl,iuby,1) = ipar0
      MolTemplateUByPar(itmpl,iuby,2) = ipar1

      return
      end
