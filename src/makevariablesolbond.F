      subroutine makeVariableSolBond(isol,zat1,zat2,idelta)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'enums.inc'
c#############################################################
c
c     This routine makes the bond between atoms zat1 and zat2
c     in solute isol flexible (if it is defined in the template
c     zmatrix
c
c     (C) Christopher Woods    19 February 2002
c
c#############################################################

      integer isol,zat1,zat2
      double precision idelta
      character*10 inam1,inam2,getSolAtomName
      
c     return if equal!
      if (zat1.eq.zat2) return
      
c     sort the bond
      call sortZBond(zat1,zat2)
      
c     now check that this is a zmatrix defined bond...
      if (ZSolMat(isol,zat1,1).eq.zat2) then
        ZSolMov(isol,zat1,1) = .true.
        ZSolDelta(isol,zat1,1) = idelta
      else
        inam1 = getSolAtomName(isol,zat1)
        inam2 = getSolAtomName(isol,zat2)
        
        write(printstring,*) "Trying to make a non-template defined internal",
     .       " bond flexible! ",inam1," - ",inam2," in solute ",isol
        call printLine(WARNING,printstring)
      endif
      return
      end
