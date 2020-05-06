      subroutine addProteinDihedral(iprot,iat1,iat2,iat3,iat4,dummy,flex,
     .                              delta,ipar0,ipar1)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'proteins.inc'
c##################################################################
c
c     This routine adds the dihedral angle between 
c     zat1-zat2-zat3-zat4 in protein iprot
c
c     (C) Christopher Woods       30 January 2002
c
c##################################################################

      integer zat1,zat2,zat3,zat4,iprot
      integer iat1,iat2,iat3,iat4
      integer iparam,getDihedral,setDirectDihedral
      integer ires1,ires2,ires3,ires4
      logical flex,dummy
      integer ipar0,ipar1
      double precision delta
      logical hasDihedral
      
      zat1 = iat1
      zat2 = iat2
      zat3 = iat3
      zat4 = iat4
      
      if (hasDihedral(iprot,zat1,zat2,zat3,zat4)) then
        return
      endif
      
c     *do not sort* if the first atom is the O, and the last is 
c     the N, since these are defined in the 'wrong' order...
      if (Zid2ResID(iprot,zat1,2).ne.-4  .and. Zid2ResID(iprot,zat4,2).ne.-1) then
        call sortZDihedral(zat1,zat2,zat3,zat4)
      endif

      ires1 = Zid2ResID(iprot,zat1,1)
      ires2 = Zid2ResID(iprot,zat2,1)
      ires3 = Zid2ResID(iprot,zat3,1)
      ires4 = Zid2ResID(iprot,zat4,1)
      
c     now add the angle...
      NDihedrals(iprot) = NDihedrals(iprot) + 1
      if (NDihedrals(iprot).gt.MAXDIHEDRALS) then
        write(printstring,*) "Exceeded maximum number of protein dihedrals! (",MAXDIHEDRALS,")"
        call printLine(FATAL,printstring)
        write(printstring,*)"Increase MAXDIHEDRALS and recompile"
        call closeProgram(-1,printstring)
        NDihedrals(iprot) = NDihedrals(iprot) - 1
        return
      endif
      Dihedral(iprot,NDihedrals(iprot),1) = zat1
      Dihedral(iprot,NDihedrals(iprot),2) = zat2
      Dihedral(iprot,NDihedrals(iprot),3) = zat3
      Dihedral(iprot,NDihedrals(iprot),4) = zat4
      
      DihedralDummy(iprot,NDihedrals(iprot)) = .false.
      if (flex) then
        FlexDihedral(iprot,NDihedrals(iprot)) = .true.
        call makeVariableDihedral(iprot,zat1,zat2,zat3,zat4,delta)
      endif

      if (dummy) then
        iparam = 0
        DihedralDummy(iprot,NDihedrals(iprot)) = .true.
      else if (ipar0.ge.0 .and. ipar1.ge.0) then
        iparam = setDirectDihedral(ipar0,ipar1)
      else
        iparam = getDihedral(ZPar(iprot,zat1),ZPar(iprot,zat2),
     .                       ZPar(iprot,zat3),ZPar(iprot,zat4))
      endif
      Dihedral(iprot,NDihedrals(iprot),5) = iparam
      
c     now add the Dihedral to the residue
      NResDihedrals(iprot,ires1) = NResDihedrals(iprot,ires1) + 1
      if (NResDihedrals(iprot,ires1).gt.MAXDIHEDRALSPERRESIDUE) then
        write(printstring,*) "Too many dihedrals in residue ",ResNam(iprot,ires1),"-",ResNum(iprot,ires1)
        call printLine(FATAL,printstring)
        write(printstring,*)"Increase MAXDIHEDRALSPERRESIDUE and recompile"
        call closeProgram(-1,printstring)
        NResDihedrals(iprot,ires1) = NResDihedrals(iprot,ires1) - 1
        return
      endif
      ResidueDihedral(iprot,ires1,NResDihedrals(iprot,ires1)) = NDihedrals(iprot)
      if (ires2.ne.ires1) then
        NResDihedrals(iprot,ires2) = NResDihedrals(iprot,ires2) + 1
        if (NResDihedrals(iprot,ires2).gt.MAXDIHEDRALSPERRESIDUE) then
          write(printstring,*) "Too many dihedrals in residue ",ResNam(iprot,ires2),"-",ResNum(iprot,ires2)
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXDIHEDRALSPERRESIDUE and recompile"
          call closeProgram(-1,printstring)
          NResDihedrals(iprot,ires2) = NResDihedrals(iprot,ires2) - 1
          return
        endif
        ResidueDihedral(iprot,ires2,NResDihedrals(iprot,ires2)) = NDihedrals(iprot)
      endif
      if (ires3.ne.ires2 .and. ires3.ne.ires1) then
        NResDihedrals(iprot,ires3) = NResDihedrals(iprot,ires3) + 1
        if (NResDihedrals(iprot,ires3).gt.MAXDIHEDRALSPERRESIDUE) then
          write(printstring,*) "Too many dihedrals in residue ",ResNam(iprot,ires3),"-",ResNum(iprot,ires3)
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXDIHEDRALSPERRESIDUE and recompile"
          call closeProgram(-1,printstring)
          NResDihedrals(iprot,ires3) = NResDihedrals(iprot,ires3) - 1
          return
        endif
        ResidueDihedral(iprot,ires3,NResDihedrals(iprot,ires3)) = NDihedrals(iprot)
      endif
      if (ires4.ne.ires3 .and. ires4.ne.ires2 .and. ires4.ne.ires1) then
        NResDihedrals(iprot,ires4) = NResDihedrals(iprot,ires4) + 1
        if (NResDihedrals(iprot,ires4).gt.MAXDIHEDRALSPERRESIDUE) then
          write(printstring,*) "Too many dihedrals in residue ",ResNam(iprot,ires4),"-",ResNum(iprot,ires4)
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXDIHEDRALSPERRESIDUE and recompile"
          call closeProgram(-1,printstring)
          NResDihedrals(iprot,ires4) = NResDihedrals(iprot,ires4) - 1
          return
        endif
        ResidueDihedral(iprot,ires4,NResDihedrals(iprot,ires4)) = NDihedrals(iprot)
      endif

c     now record the residues involved in the dihedral
      Dihedral(iprot,NDihedrals(iprot),6) = ires1
      Dihedral(iprot,NDihedrals(iprot),7) = ires2
      Dihedral(iprot,NDihedrals(iprot),8) = ires3
      Dihedral(iprot,NDihedrals(iprot),9) = ires4

      return
      
      end
