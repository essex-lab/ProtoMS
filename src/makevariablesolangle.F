      subroutine makeVariableSolAngle(isol,zat1,zat2,zat3,idelta)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'enums.inc'
c#######################################################
c
c     This makes the template defined solute angle flexible
c
c     (C) Christopher Woods   19 February 2002
c
c#######################################################

      integer isol,zat1,zat2,zat3
      double precision idelta
      character*10 inam1,inam2,inam3,getSolAtomName
      
      if (zat1.eq.zat3) return
      call sortZAngle(zat1,zat2,zat3)
      
      if (ZSolMat(isol,zat1,1).eq.zat2 .and. ZSolMat(isol,zat1,2).eq.zat3) then
        ZSolMov(isol,zat1,2) = .true.
        ZSolDelta(isol,zat1,2) = idelta
      else
        inam1 = getSolAtomName(isol,zat1)
        inam2 = getSolAtomName(isol,zat2)
        inam3 = getSolAtomName(isol,zat3)
                
        write(printstring,*) "Trying to make a non-template defined internal",
     .       " angle flexible! ",inam1," - ",inam2," - ",inam3," in solute ",isol
        call printLine(WARNING,printstring)
      endif
      return
      end
