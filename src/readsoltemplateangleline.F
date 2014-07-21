      subroutine readSolTemplateAngleLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c###########################################################
c
c     This routine reads in the angles for a molecule template
c
c     (C) Christopher Woods
c
c###########################################################

      character*(*) line
      integer itmpl
      character*4 iatm1,ires1,iatm2,ires2,iatm3,ires3
      logical dummy,flex
      double precision delta
      integer ipar0,ipar1
      integer findMolTemplateAngle
      integer iang
      double precision deg2Rad

c     format is 'angle iatm1 ires1 iatm2 ires2 iatm3 ires3 (dummy) (flex delta) (param ipar0 ipar1)

      if (nwords.lt.7) then
        call printLine(WARNING,"Incorrectly formatted solute template angle line?")
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
      
c     see if the angle already exists      
      iang = findMolTemplateAngle(iatm1,ires1,iatm2,ires2,iatm3,ires3,itmpl)
      if (iang.ne.0) then
        write(printstring,*) "Overwriting information for solute angle (",
     .     iatm1,",",ires1,")-(",iatm2,",",ires2,")-(",iatm3,",",ires3,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
      else
        iang = NMolTemplateAngles(itmpl) + 1
        if (iang.gt.MAXMOLTEMPLATEANGLES) then
          write(printstring,*) "Exceeded maximum number of molecule template angles (",
     .                          MAXMOLTEMPLATEANGLES,")"
          call printLine(WARNING,printstring)
          call printLine(WARNING,"Increase MAXMOLTEMPLATEANGLES and recompile")
          return
        endif
        NMolTemplateAngles(itmpl) = iang
      endif
      
      MolTemplateAngName(itmpl,iang,1) = iatm1
      MolTemplateAngRes(itmpl,iang,1) = ires1
      MolTemplateAngName(itmpl,iang,2) = iatm2
      MolTemplateAngRes(itmpl,iang,2) = ires2
      MolTemplateAngName(itmpl,iang,3) = iatm3
      MolTemplateAngRes(itmpl,iang,3) = ires3
      MolTemplateAngleType(itmpl,iang,1)= dummy
      MolTemplateAngleType(itmpl,iang,2)= flex
      MolTemplateAngleDelta(itmpl,iang) = deg2Rad(delta)
      MolTemplateAngPar(itmpl,iang,1) = ipar0
      MolTemplateAngPar(itmpl,iang,2) = ipar1

      return
      end
