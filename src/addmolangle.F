      subroutine addMolAngle(isol,tmol,iang)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'proteins.inc'
      include 'enums.inc'
c###########################################################
c
c     Angle equivelent of addMolBond
c
c     (C) Christopher Woods    18 February 2002
c
c###########################################################

      integer isol,tmol,iang
      integer zat1,zat2,zat3
      integer findSoluteAtom
      
      zat1 = findSoluteAtom(isol,MolTemplateAngName(tmol,iang,1),MolTemplateAngRes(tmol,iang,1))
      zat2 = findSoluteAtom(isol,MolTemplateAngName(tmol,iang,2),MolTemplateAngRes(tmol,iang,2))
      zat3 = findSoluteAtom(isol,MolTemplateAngName(tmol,iang,3),MolTemplateAngRes(tmol,iang,3))
      
      if (zat1.eq.0 .or. zat2.eq.0 .or. zat3.eq.0) then
        write(printstring,*) "Could not find solute angle in solute ",isol," ",
     .      MolTemplateAngName(tmol,iang,1),MolTemplateAngRes(tmol,iang,1)," ",
     .      MolTemplateAngName(tmol,iang,2),MolTemplateAngRes(tmol,iang,2)," ",
     .      MolTemplateAngName(tmol,iang,3),MolTemplateAngRes(tmol,iang,3)," ",
     .      zat1,zat2,zat3
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Check the template file to see if the atoms really exist")
        return
      endif
      
      call addSoluteAngle(isol,zat1,zat2,zat3,MolTemplateAngleType(tmol,iang,1),
     .  MolTemplateAngleType(tmol,iang,2),MolTemplateAngleDelta(tmol,iang),
     .  MolTemplateAngPar(tmol,iang,1),MolTemplateAngPar(tmol,iang,2))

      return
      end
