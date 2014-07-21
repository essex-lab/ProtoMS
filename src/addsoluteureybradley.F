      subroutine addSoluteUreyBradley(isol,zat1,zat2,zat3,ipar0,ipar1)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'solutes.inc'
c####################################################################
c
c     This routine adds a Urey-Bradley term in solute isol, between atoms with
c     zmatrix id's zat1,zat2,zat3
c
c     (C) Christopher Woods
c
c####################################################################

      integer isol,zat1,zat2,zat3,nubys,iparam
      integer getUreyBradley
      integer setDirectUreyBradley
      integer ipar0,ipar1
      logical hasSoluteUreyBradley

c     make sure that the Urey-Bradley term has not already been added
      if (hasSoluteUreyBradley(isol,zat1,zat2,zat3)) return

c     now add the uby
      NSolUreyBradleys(isol) = NSolUreyBradleys(isol) + 1
      nubys = NSolUreyBradleys(isol)
      if (nubys.gt.MAXSOLUTEUREYBRADLEYS) then
        write(printstring,*) "Exceeded maximum number of solute Urey-Bradley terms! (",
     .                        MAXSOLUTEUREYBRADLEYS,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Increase MAXSOLUTEUREYBRADLEYS and recompile")
        NSolUreyBradleys(isol) = NSolUreyBradleys(isol) - 1
        return
      endif
      
      SolUreyBradley(isol,nubys,1) = zat1
      SolUreyBradley(isol,nubys,2) = zat2
      SolUreyBradley(isol,nubys,3) = zat3
      
      if (zat1.le.0 .or. zat2.le.0 .or. zat3.le.0) then
c       this uby involves one of the three special dummy atoms - this must thus
c       be a dummy uby
        iparam = 0
      else if (ipar0.ge.0 .and. ipar1.ge.0) then
        iparam = setDirectUreyBradley(ipar0,ipar1)
      else
        iparam = getUreyBradley(ZSolPar(isol,zat1),ZSolPar(isol,zat2),
     .                          ZSolPar(isol,zat3))
      endif
      
      SolUreyBradley(isol,nubys,4) = iparam
      return
      end
