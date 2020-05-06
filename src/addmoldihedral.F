      subroutine addMolDihedral(isol,tmol,idih)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'proteins.inc'
      include 'enums.inc'
c############################################################
c
c     Dihedral equivelent of addMolBond
c
c     Christopher Woods
c
c############################################################

      integer isol,tmol,idih,zat1,zat2,zat3,zat4
      integer findSoluteAtom

      zat1 = findSoluteAtom(isol,MolTemplateDihName(tmol,idih,1),MolTemplateDihRes(tmol,idih,1))
      zat2 = findSoluteAtom(isol,MolTemplateDihName(tmol,idih,2),MolTemplateDihRes(tmol,idih,2))
      zat3 = findSoluteAtom(isol,MolTemplateDihName(tmol,idih,3),MolTemplateDihRes(tmol,idih,3))
      zat4 = findSoluteAtom(isol,MolTemplateDihName(tmol,idih,4),MolTemplateDihRes(tmol,idih,4))
      
      if (zat1.eq.0 .or. zat2.eq.0 .or. zat3.eq.0 .or. zat4.eq.0) then
        write(printstring,*) "Could not find dihedral in solute ",isol," ",
     .      MolTemplateDihName(tmol,idih,1),MolTemplateDihRes(tmol,idih,1)," ",
     .      MolTemplateDihName(tmol,idih,2),MolTemplateDihRes(tmol,idih,2)," ",
     .      MolTemplateDihName(tmol,idih,3),MolTemplateDihRes(tmol,idih,3)," ",
     .      MolTemplateDihName(tmol,idih,4),MolTemplateDihRes(tmol,idih,4)," ",
     .      zat1,zat2,zat3,zat4
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Check the template file to see if the atoms really exist")
        return
      endif
      
      call addSoluteDihedral(isol,zat1,zat2,zat3,zat4,
     .  MolTemplateDihedralType(tmol,idih,1),MolTemplateDihedralType(tmol,idih,2),
     .  MolTemplateDihedralDelta(tmol,idih),
     .  MolTemplateDihPar(tmol,idih,1),MolTemplateDihPar(tmol,idih,2),
     .  MolTemplateDihNWells(tmol,idih),
     .  MolTemplateDihWell(tmol,idih,1),MolTemplateDihWell(tmol,idih,2),
     .  MolTemplateDihWell(tmol,idih,3),MolTemplateDihWell(tmol,idih,4))

      return
      end
