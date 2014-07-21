
      subroutine gbsavesystem
      include 'dimensions.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'printer.inc'
      include 'simparams.inc'
      include 'enums.inc'
      include 'flags.inc'      
      include 'gb.inc'
c#####################################
c     Save the coordinates and centers 
c     of the moving bits of a move
c
c     Julien Michel November 2004
c#####################################
      integer i,j,nres
      logical hasSoluteMoved,hasResidueMoved,hasSolventMoved

      if (NSolutes.gt. ZERO) then
        do 20 i=1,NSolutes
          MovedSolute(i) = 0
          if (hasSoluteMoved(i)) then
c            write(printstring,*) 'SOLUTE ',i,' has moved'
c            call printLine(DEBUG,printstring)
c           Add to list of moved solutes
            call SaveSoluteGB(i)
            MovedSolute(i) = countGBSolute          
          endif
 20     continue
      endif

      if (NProteins.gt. ZERO) then
        do 30 i=1,NProteins
          nres = NResidues(i)
          do 40 j=1,nres
          MovedResidue(i,j) = 0
          if (hasResidueMoved(i,j)) then
c           Add to list of moved solutes
            call SaveResidueGB(i,j)
            MovedResidue(i,j) = countGBResidue          
          endif
 40       continue
 30     continue
      endif
      
      if (NSolvents .gt. ZERO) then
        do 50 i=1,NSolvents
          MovedSolvent(i) = 0
          if (hasSolventMoved(i)) then
c           Add to list of moved solutes
            call SaveSolventGB(i)
            MovedSolvent(i) = countGBSolvent
          endif
 50     continue
      endif

      return
      end
