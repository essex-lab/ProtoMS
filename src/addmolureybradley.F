      subroutine addMolUreyBradley(isol,tmol,iuby)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'enums.inc'
c#########################################################
c
c     This routine adds the template urey-bradley iuby, from the 
c     molecule template tmol, to the solute isol
c
c     (C) Christopher Woods
c
c##########################################################

      integer isol,tmol,iuby,zat1,zat2,zat3
      integer findSoluteAtom
      
      zat1 = findSoluteAtom(isol,MolTemplateUByName(tmol,iuby,1),MolTemplateUByRes(tmol,iuby,1))
      zat2 = findSoluteAtom(isol,MolTemplateUByName(tmol,iuby,2),MolTemplateUByRes(tmol,iuby,2))
      zat3 = findSoluteAtom(isol,MolTemplateUByName(tmol,iuby,3),MolTemplateUByRes(tmol,iuby,3))
      if (zat1.eq.0 .or. zat2.eq.0 .or. zat3.eq.0) then
        write(printstring,*) "Could not find Urey-Bradley term in solute ",isol," ",
     .                MolTemplateUByName(tmol,iuby,1),MolTemplateUByRes(tmol,iuby,1)," ",
     .                MolTemplateUByName(tmol,iuby,2),MolTemplateUByRes(tmol,iuby,2)," ",
     .                zat1,zat2
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Check the template file to see if the atoms really exist")
        return
      endif
      
      call addSoluteUreyBradley(isol,zat1,zat2,zat3,
     .         MolTemplateUByPar(tmol,iuby,1),MolTemplateUByPar(tmol,iuby,2))
      return
      end
