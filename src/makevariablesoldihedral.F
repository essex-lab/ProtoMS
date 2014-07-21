      subroutine makeVariableSolDihedral(isol,zat1,zat2,zat3,zat4,idelta)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'enums.inc'
c##############################################################
c
c     This makes this solute dihedral flexible (if defined
c     in the template file)
c
c     (C) Christopher Woods    19 February 2002
c
c##############################################################

      integer isol,zat1,zat2,zat3,zat4
      double precision idelta
      character*10 inam1,inam2,inam3,inam4,getSolAtomName
      
      call sortZDihedral(zat1,zat2,zat3,zat4)
      
      if (ZSolMat(isol,zat1,1).eq.zat2 .and. ZSolMat(isol,zat1,2).eq.zat3
     .                    .and. ZSolMat(isol,zat1,3).eq.zat4) then
        ZSolMov(isol,zat1,3) = .true.
        ZSolDelta(isol,zat1,3) = idelta
      else
        inam1 = getSolAtomName(isol,zat1)
        inam2 = getSolAtomName(isol,zat2)
        inam3 = getSolAtomName(isol,zat3)
        inam4 = getSolAtomName(isol,zat4)        
        
        write(printstring,*) "Trying to make a non-template defined internal",
     .       " dihedral flexible! ",inam1," - ",inam2," - ",inam3," - ",inam4," in solute ",isol
        call printLine(WARNING,printstring)
      endif


      return
      end
