      subroutine addSoluteDihedral(isol,zat1,zat2,zat3,zat4,dummy,flex,delta,
     .                                ipar0,ipar1,nwells,well1,well2,well3,well4)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'solutes.inc'
c###################################################################
c
c     This routine adds the solute dihedral angle to solute isol,
c     between atoms zat1,zat2,zat3,zat4, with type ityp, flex iflx
c     and delta idelta
c
c     (C) Christopher Woods    19 February 2002
c
c###################################################################

      integer isol,zat1,zat2,zat3,zat4,ndihs,iparam
      integer getDihedral
      integer setDirectDihedral
      integer ipar0,ipar1
      logical flex,dummy
      double precision delta
      logical hasSoluteDihedral
      integer nwells
      double precision well1,well2,well3,well4

      if (hasSoluteDihedral(isol,zat1,zat2,zat3,zat4)) return
      
      NSolDihedrals(isol) = NSolDihedrals(isol) + 1
      ndihs = NSolDihedrals(isol)
      if (ndihs.gt.MAXSOLUTEDIHEDRALS) then
        write(printstring,*) "Exceeded maximum number of solute dihedrals! (",
     .                        MAXSOLUTEDIHEDRALS,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Increase MAXSOLUTEDIHEDRALS and recompile")
        NSolDihedrals(isol) = NSolDihedrals(isol) - 1
        return
      endif
      
      SolDihedral(isol,ndihs,1) = zat1
      SolDihedral(isol,ndihs,2) = zat2
      SolDihedral(isol,ndihs,3) = zat3
      SolDihedral(isol,ndihs,4) = zat4
      
      if (flex) then
        FlexSolDihedral(isol,ndihs) = .true.
        call makeVariableSolDihedral(isol,zat1,zat2,zat3,zat4,delta)
        if (nwells .gt. 0) then
          ZSolDihNWells(isol,zat1) = nwells
          ZSolDihWell(isol,zat1,1) = well1
          ZSolDihWell(isol,zat1,2) = well2
          ZSolDihWell(isol,zat1,3) = well3
          ZSolDihWell(isol,zat1,4) = well4
c          write(printstring,*) ' Well on sol ',isol,' at ',zat1,
c     .      ' nwells ',nwells,' well2 ',ZSolDihWell(isol,zat1,2)
c          call printLine(DEBUG,printstring)
        endif
      endif
      
      if (dummy .or. zat1.le.0 .or. zat2.le.0 .or. zat3.le.0 .or. zat4.le.0) then
        iparam = 0
        SolDihedralDummy(isol,ndihs) = .true.
      else if (ipar0.ge.0 .and. ipar1.ge.0) then
        iparam = setDirectDihedral(ipar0,ipar1)
      else
        iparam = getDihedral(ZSolPar(isol,zat1),ZSolPar(isol,zat2),
     .                       ZSolPar(isol,zat3),ZSolPar(isol,zat4))
      endif
      
      SolDihedral(isol,ndihs,5) = iparam


      return
      end
