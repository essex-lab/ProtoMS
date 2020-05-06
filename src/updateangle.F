      subroutine updateAngle(iprot,ires,iresang)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'enums.inc'
c###########################################################################
c
c     This routine updates the angle parameters for angle iresang in 
c     residue ires, protein iprot.
c
c     (C) Christopher Woods  4 February 2002
c
c###########################################################################

      integer iprot,ires,iresang,iang
      integer zat1,zat2,zat3
      integer findAngleParam
      
      if (iresang.le.0 .or. iresang.gt.NResAngles(iprot,ires)) return
      
c     get the angle number in the protein of this bond...
      iang = ResidueAngle(iprot,ires,iresang)
c     if this is a dummy angle, then return
      if (Angle(iprot,iang,4).lt.0) then
        write(printstring,*) "Attempting to convert residue with variable angle parameters! - ",iprot,ires
        call closeProgram(-1,printstring)
        return
      endif
      if (AngleDummy(iprot,iang)) return
      zat1 = Angle(iprot,iang,1)
      zat2 = Angle(iprot,iang,2)
      zat3 = Angle(iprot,iang,3)

      Angle(iprot,iang,4) = findAngleParam(ZPar(iprot,zat1),ZPar(iprot,zat2),ZPar(iprot,zat3))
      
      return
      
      end
