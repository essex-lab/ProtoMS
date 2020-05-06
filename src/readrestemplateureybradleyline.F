      subroutine readResTemplateUreyBradleyLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c###########################################################
c
c     This routine reads in a residue Urey-Bradley template line
c
c     (C) Christopher Woods
c
c###########################################################

      character*(*) line
      integer itmpl
      character*4 at1,at2,at3
      logical flex,dummy
      integer ipar0,ipar1,iuby,findTemplateUreyBradley
      double precision delta
      
c     format is "ureybradley atm1 atm2 atm3 (dummy) (param ipar0 ipar1)
c      (ubys cannot be flexible)

      if (nwords.lt.4) then
        call printLine(WARNING,"Incorrectly formatted residue Urey-Bradley template line?")
        call printLine(WARNING,line)
        return
      endif
      
c     read in the atoms involved
      at1 = words(2)(1:4)
      at2 = words(3)(1:4)
      at3 = words(4)(1:4)

      call sortTemplateUreyBradley(at1,at2,at3)
      iuby = findTemplateUreyBradley(at1,at2,at3,itmpl)
      if (iuby.eq.0) then
        iuby = NTemplateUreyBradleys(itmpl) + 1
        if (iuby.gt.MAXTEMPLATEUREYBRADLEYS) then
          write(printstring,*) "Exceeded maximum number of residue template Urey-Bradley terms (",
     .                         MAXTEMPLATEUREYBRADLEYS,")"
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXTEMPLATEUREYBRADLEYS and recompile"
          call closeProgram(-1,printstring)
          return
        endif
        NTemplateUreyBradleys(itmpl) = iuby
      else
        write(printstring,*) "Overwriting residue template Urey-Bradley for ",at1,"-",
     .                       at2,"-",at3
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
      endif
            
      TemplateUByName(itmpl,iuby,1) = at1
      TemplateUByName(itmpl,iuby,2) = at2
      TemplateUByName(itmpl,iuby,3) = at3
      
      call getDOFInfo(dummy,flex,delta,ipar0,ipar1)
      TemplateUreyBradleyType(itmpl,iuby) = dummy
      if (flex) then
        call printLine(WARNING,"Cannot have flexible Urey-Bradleys!")
        call printLine(WARNING,line)
      endif
      
      TemplateUreyBradleyPar(itmpl,iuby,1) = ipar0
      TemplateUreyBradleyPar(itmpl,iuby,2) = ipar1
      
      return
      end  
