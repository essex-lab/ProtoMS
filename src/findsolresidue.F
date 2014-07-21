      integer function findSolResidue(isol,iresnum,iresnam)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'solutes.inc'
      include 'enums.inc'
c#############################################################
c
c     This function sees if solute isol already contains
c     a residue numbered 'iresnum' or called 'iresnam', 
c     and returns the index to the residue if it does.
c
c     It returns 0 if the residue cannot be found
c
c     (C) Christopher Woods, October 2004
c
c#############################################################

      integer isol,iresnum,i
      character*4 iresnam
      
      findSolResidue = 0
      
      if (isol.le.0 .or. isol.gt.NSolutes) return
      
      if (NSolResidues(isol).le.0) return
      
      do i=1,NSolResidues(isol)
        if (SolResNam(isol,i).eq.iresnam) then
c          if (SolResNum(isol,i).ne.iresnum) call printLine(WARNING,"Serious problem with solute pdb!")
          findSolResidue = i
          return
        else if (SolResNum(isol,i).eq.iresnum) then
c          if (SolResNam(isol,i).ne.iresnam) call printLine(WARNING,"Serious problem with solute pdb!")
          findSolResidue = i
          return
        endif
      enddo
      
      return
      end
