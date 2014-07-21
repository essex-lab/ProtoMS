      subroutine addProteinAngle(iprot,iat1,iat2,iat3,dummy,flex,delta,ipar0,ipar1)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'proteins.inc'
c################################################################
c
c     This adds the angle in iprot between zat1-zat2-zat3
c
c     (C) Christopher Woods         30 January 2002
c
c###############################################################

      integer iprot,zat1,zat2,zat3
      integer iat1,iat2,iat3
      integer ires1,ires2,ires3,iparam
      integer getAngle,setDirectAngle
      logical flex,dummy
      integer ipar0,ipar1
      double precision delta
      logical hasAngle
      
      zat1 = iat1
      zat2 = iat2
      zat3 = iat3
      
      if (hasAngle(iprot,zat1,zat2,zat3)) then
        return
      endif
      
      call sortZAngle(zat1,zat2,zat3)
      
      ires1 = Zid2ResID(iprot,zat1,1)
      ires2 = Zid2ResID(iprot,zat2,1)
      ires3 = Zid2ResID(iprot,zat3,1)
      if (ires1.le.0 .or. ires2.le.0 .or. ires3.le.0) then
        write(printstring,10) ires1,ires2,ires3
10      format("There is a problem with Zid2ResID, ",i3,", ",i3," and ",i3)
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is a probably a program bug or template typo?")
        return
      endif

c     now add the angle...
      NAngles(iprot) = NAngles(iprot) + 1
      if (NAngles(iprot).gt.MAXANGLES) then
        write(printstring,*) "Exceeded maximum number of protein angles! (",MAXANGLES,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Increase MAXANGLES and recompile")
        NAngles(iprot) = NAngles(iprot) - 1
        return
      endif
      Angle(iprot,NAngles(iprot),1) = zat1
      Angle(iprot,NAngles(iprot),2) = zat2
      Angle(iprot,NAngles(iprot),3) = zat3
      AngleDummy(iprot,NAngles(iprot)) = .false.
      if (flex) then
        FlexAngle(iprot,NAngles(iprot)) = .true.
        call makeVariableAngle(iprot,zat1,zat2,zat3,delta)
      endif
      
      if (dummy) then
c       this is a dummy angle
        iparam = 0
        AngleDummy(iprot,NAngles(iprot)) = .true.
      else if (ipar0.ge.0 .and. ipar1.ge.0) then
        iparam = setDirectAngle(ipar0,ipar1)
      else
         iparam = getAngle(ZPar(iprot,zat1),ZPar(iprot,zat2),ZPar(iprot,zat3))
      endif
      Angle(iprot,NAngles(iprot),4) = iparam

c     now add the Angle to the residue
      NResAngles(iprot,ires1) = NResAngles(iprot,ires1) + 1
      if (NResAngles(iprot,ires1).gt.MAXANGLESPERRESIDUE) then
        write(printstring,*) "Too many angles in residue ",ires1
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Increase MAXANGLESPERRESIDUE and recompile")
        NResAngles(iprot,ires1) = NResAngles(iprot,ires1) - 1
        return
      endif
      ResidueAngle(iprot,ires1,NResAngles(iprot,ires1)) = NAngles(iprot)
      if (ires2.ne.ires1) then
        NResAngles(iprot,ires2) = NResAngles(iprot,ires2) + 1
        if (NResAngles(iprot,ires2).gt.MAXANGLESPERRESIDUE) then
          write(printstring,*) "Too many angles in residue ",ires2
          call printLine(WARNING,printstring)
          call printLine(WARNING,"Increase MAXANGLESPERRESIDUE and recompile")
          NResAngles(iprot,ires2) = NResAngles(iprot,ires2) - 1
          return
        endif
        ResidueAngle(iprot,ires2,NResAngles(iprot,ires2)) = NAngles(iprot)
      endif
      if (ires3.ne.ires1 .and. ires3.ne.ires2) then
        NResAngles(iprot,ires3) = NResAngles(iprot,ires3) + 1
        if (NResAngles(iprot,ires3).gt.MAXANGLESPERRESIDUE) then
          write(printstring,*) "Too many angles in residue ",ires3
          call printLine(WARNING,printstring)
          call printLine(WARNING,"Increase MAXANGLESPERRESIDUE and recompile")
          NResAngles(iprot,ires3) = NResAngles(iprot,ires3) - 1
          return
        endif
        ResidueAngle(iprot,ires3,NResAngles(iprot,ires3)) = NAngles(iprot)
      endif
      
c     now record the residues involved in the angle
      Angle(iprot,NAngles(iprot),5) = ires1
      Angle(iprot,NAngles(iprot),6) = ires2
      Angle(iprot,NAngles(iprot),7) = ires3
      
      return
      
      end
