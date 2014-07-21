      subroutine addSoluteBond(isol,zat1,zat2,dummy,flex,delta,ipar0,ipar1)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'solutes.inc'
c####################################################################
c
c     This routine adds a bond in solute isol, between atoms with
c     zmatrix id's zat1,zat2, with type ityp, flex iflx and
c     delta idelta
c
c     (C) Christopher Woods  18 February 2002
c
c####################################################################

      integer isol,zat1,zat2,nbnds,iparam
      integer getBond
      integer setDirectBond
      logical flex,dummy
      integer ipar0,ipar1
      double precision delta
      logical hasSoluteBond

c     make sure that the bond has not yet been added
      if (hasSoluteBond(isol,zat1,zat2)) return

c     now add the bond
      NSolBonds(isol) = NSolBonds(isol) + 1
      nbnds = NSolBonds(isol)
      if (nbnds.gt.MAXSOLUTEBONDS) then
        write(printstring,*) "Exceeded maximum number of solute bonds! (",MAXSOLUTEBONDS,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Increase MAXSOLUTEBONDS and recompile")
        NSolBonds(isol) = NSolBonds(isol) - 1
        return
      endif
      
      SolBond(isol,nbnds,1) = zat1
      SolBond(isol,nbnds,2) = zat2

      if (flex) then
        FlexSolBond(isol,nbnds) = .true.
        call makeVariableSolBond(isol,zat1,zat2,delta)
      endif

c     if this has been explicitly set as a dummy, or is a bond back to the
c     three original dummies, then this is a dummy bond
      if (dummy .or. zat1.le.0 .or. zat2.le.0) then
        iparam = 0
        SolBondDummy(isol,nbnds) = .true.
      else if (ipar0.ge.0 .and. ipar1.ge.0) then
        iparam = setDirectBond(ipar0,ipar1)
      else
        iparam = getBond(ZSolPar(isol,zat1),ZSolPar(isol,zat2))
      endif
      
      SolBond(isol,nbnds,3) = iparam
      return
      end
