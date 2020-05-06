      subroutine addMolBond(isol,tmol,ibnd)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'proteins.inc'
      include 'enums.inc'
c#########################################################
c
c     This routine adds the template bond ibnd, from the 
c     molecule template tmol, to the solute isol
c
c     (C) Christopher Woods   18 February 2002
c
c##########################################################

      integer isol,tmol,ibnd,zat1,zat2
      integer findSoluteAtom
      
      zat1 = findSoluteAtom(isol,MolTemplateBndName(tmol,ibnd,1),MolTemplateBndRes(tmol,ibnd,1))
      zat2 = findSoluteAtom(isol,MolTemplateBndName(tmol,ibnd,2),MolTemplateBndRes(tmol,ibnd,2))
      if (zat1.eq.0 .or. zat2.eq.0) then
        write(printstring,*) "Could not find solute bond in solute ",isol," ",
     .                MolTemplateBndName(tmol,ibnd,1),MolTemplateBndRes(tmol,ibnd,1)," ",
     .                MolTemplateBndName(tmol,ibnd,2),MolTemplateBndRes(tmol,ibnd,2)," ",
     .                zat1,zat2
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Check the template file to see if the atoms really exist")
        return
      endif
      
      call addSoluteBond(isol,zat1,zat2,MolTemplateBondType(tmol,ibnd,1),
     .         MolTemplateBondType(tmol,ibnd,2),MolTemplateBondDelta(tmol,ibnd),
     .         MolTemplateBndPar(tmol,ibnd,1),MolTemplateBndPar(tmol,ibnd,2))
      return
      end
