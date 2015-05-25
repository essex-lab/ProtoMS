      subroutine readResTemplateAngleLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c##########################################################
c
c     This routine reads in a residue angle template line
c
c     (C) Christopher Woods
c
c##########################################################

      character*(*) line
      integer itmpl
      character*4 at1,at2,at3
      logical flex,dummy
      integer ipar0,ipar1,iang,findTemplateAngle
      double precision delta,deg2Rad
      
c     format is "angle atm1 atm2 atm3 (flex idelta) (dummy) (param ipar0 ipar1)

      if (nwords.lt.4) then
        call printLine(WARNING,"Incorrectly formatted residue angle template line?")
        call printLine(WARNING,line)
        return
      endif
      
c     read in the atoms involved
      at1 = words(2)(1:4)
      at2 = words(3)(1:4)
      at3 = words(4)(1:4)
      
      call sortTemplateAngle(at1,at2,at3)
      
      iang = findTemplateAngle(at1,at2,at3,itmpl)
      if (iang.eq.0) then
        iang = NTemplateAngles(itmpl) + 1
        if (iang.gt.MAXTEMPLATEANGLES) then
          write(printstring,*) "Exceeded maximum number of residue template angles (",
     .                          MAXTEMPLATEANGLES,")"
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXTEMPLATEANGLES and recompile"
          call closeProgram(-1,printstring)
          return
        endif
        NTemplateAngles(itmpl) = iang
      else
        write(printstring,*) "Overwriting residue template angle for ",at1,"-",at2,
     .                             "-",at3
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
      endif
      
      TemplateAngName(itmpl,iang,1) = at1
      TemplateAngName(itmpl,iang,2) = at2
      TemplateAngName(itmpl,iang,3) = at3
      
      call getDOFInfo(dummy,flex,delta,ipar0,ipar1)
      TemplateAngleType(itmpl,iang,1) = dummy
      TemplateAngleType(itmpl,iang,2) = flex
      TemplateAngleDelta(itmpl,iang) = deg2Rad(delta)
      TemplateAnglePar(itmpl,iang,1) = ipar0
      TemplateAnglePar(itmpl,iang,2) = ipar1
      
      return
      end
