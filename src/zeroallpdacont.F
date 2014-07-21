      subroutine ZeroAllPDACont
c###############################################
c     Reset ZContBorn to ZERO, during initialisation !
c
c
c     Julien Michel November 2004
c#############################################
      include 'dimensions.inc'
      include 'constants.inc'
      include 'parameters.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'solvents.inc'      
      include 'gb.inc'
      integer i,j,k,nres,kstart,kend,nats

      if (NSolutes .gt. ZERO) then
        do 10 i=1,NSolutes
          nats = NSolAtoms(i)
          do 11 j=1,nats
            ZSolContBorn(i,j) = ZERO
            ZSolContBornF(i,j) = ZERO
            ZSolContBornB(i,j) = ZERO
 11       continue
 10     continue
      endif
      
      if (NProteins .gt. ZERO) then
        do 20 i=1,NProteins
          nres = NResidues(i)
          do 21 j=1,nres
            kstart = BBZid(i,j,1)
            kend = BBZid(i,j,4)
            do 22 k=kstart,kend
              ZProtContBorn(i,k) = ZERO
              ZProtContBornF(i,k) = ZERO
              ZProtContBornB(i,k) = ZERO
 22         continue
 21       continue
 20     continue
      endif

      if (NSolvents .gt. ZERO) then
        do 30 i=1,NSolvents
          nats = NSvnAtoms(i)
          do 31 j=1,nats
            ZSvnContBorn(i,j) = ZERO
            ZSvnContBornF(i,j) = ZERO
            ZSvnContBornB(i,j) = ZERO
 31       continue
 30     continue
      endif

      return
      end
