      subroutine updateDihedral(iprot,ires,iresdih)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'flags.inc'
      include 'proteins.inc'
      include 'enums.inc'
c###########################################################################
c
c     This routine updates the dihedral parameters for dihedral iresdih in 
c     residue ires, protein iprot.
c
c     (C) Christopher Woods  4 February 2002
c
c###########################################################################

      integer iprot,ires,iresdih
      integer idih,zat1,zat2,zat3,zat4
      integer findDihedralParam
      character*13 nameString,inam1,inam2,inam3,inam4
      
      if (iresdih.le.0 .or. iresdih.gt.NResDihedrals(iprot,ires)) return
      
c     get the dihedral number in the protein of this bond...
      idih = ResidueDihedral(iprot,ires,iresdih)
      if (Dihedral(iprot,idih,5).lt.0) then
        write(printstring,*) "Attempting to convert residue with variable dihedral parameters - ",iprot,ires
        call closeProgram(-1,printstring)
        return
      endif
      if (DihedralDummy(iprot,idih)) return
      zat1 = Dihedral(iprot,idih,1)
      zat2 = Dihedral(iprot,idih,2)
      zat3 = Dihedral(iprot,idih,3)
      zat4 = Dihedral(iprot,idih,4)

      if (protodebug) then
        inam1 = nameString(iprot,zat1)
        inam2 = nameString(iprot,zat2)
        inam3 = nameString(iprot,zat3)
        inam4 = nameString(iprot,zat4)
        write(printstring,*) inam1,inam2,inam3,inam4
        call printLine(INFO,printstring)
      endif
        
      Dihedral(iprot,idih,5) = findDihedralParam(ZPar(iprot,zat1),ZPar(iprot,zat2),
     .                                           ZPar(iprot,zat3),ZPar(iprot,zat4))
      
      return
      
      end
