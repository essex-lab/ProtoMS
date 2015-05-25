      subroutine addSoluteAngle(isol,zat1,zat2,zat3,dummy,flex,delta,ipar0,ipar1)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'solutes.inc'
c##################################################################
c
c     This routine adds the solute angle to solute isol, between
c     atoms zat1,zat2,zat3, with type ityp, flex iflx, and delta idelta
c
c     (C) Christopher Woods    19 February 2002
c
c##################################################################

      integer isol,zat1,zat2,zat3,nangs,iparam
      integer getAngle
      integer setDirectAngle
      integer ipar0,ipar1
      logical dummy,flex
      double precision delta
      logical hasSoluteAngle
c      character*2 amber      

      if (hasSoluteAngle(isol,zat1,zat2,zat3)) return
      
      NSolAngles(isol) = NSolAngles(isol) + 1
      nangs = NSolAngles(isol)
      if (nangs.gt.MAXSOLUTEANGLES) then
        write(printstring,*) "Exceeded maximum number of solute angles! (",MAXSOLUTEANGLES,")"
        call printLine(FATAL,printstring)
        write(printstring,*)"Increase MAXSOLUTEANGLES and recompile"
        call closeProgram(-1,printstring)
        NSolAngles(isol) = NSolAngles(isol) - 1
        return
      endif
      
      SolAngle(isol,nangs,1) = zat1
      SolAngle(isol,nangs,2) = zat2
      SolAngle(isol,nangs,3) = zat3
      
      if (flex) then
        FlexSolAngle(isol,nangs) = .true.
        call makeVariableSolAngle(isol,zat1,zat2,zat3,delta)
      endif

      if (dummy .or. zat1.le.0 .or. zat2.le.0 .or. zat3.le.0) then
        iparam = 0
        SolAngleDummy(isol,nangs) = .true.
      else if (ipar0.ge.0 .and. ipar1.ge.0) then
        iparam = setDirectAngle(ipar0,ipar1)
      else
        iparam = getAngle(ZSolPar(isol,zat1),ZSolPar(isol,zat2),ZSolPar(isol,zat3))
      endif

      SolAngle(isol,nangs,4) = iparam

c      write(printstring,*) 'ANG BTW ',zat1,'-',amber(ZSolPar(isol,zat1)),'-',' ',zat2,'-',amber(ZSolPar(isol,zat2)),' ',zat3,'-',amber(ZSolPar(isol,zat1)),' iparam ',iparam,' ipar0 ',ipar0,' ipar1 ',ipar1,' flex ? ',flex
c      call printLine(DEBUG,printstring)
      return
      end
