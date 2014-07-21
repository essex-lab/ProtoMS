      subroutine ZeroPDACont
c###############################################
c Set the Contributions to the Effective Born radii to Zero, is called
c before computing Born radii
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
            NewZSolContBorn(i,j) = ZERO
            OldZSolContBorn(i,j) = ZERO
            NewZSolContBornF(i,j) = ZERO
            OldZSolContBornF(i,j) = ZERO
            NewZSolContBornB(i,j) = ZERO
            OldZSolContBornB(i,j) = ZERO
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
              NewZProtContBorn(i,k) = ZERO
              OldZProtContBorn(i,k) = ZERO
              NewZProtContBornF(i,k) = ZERO
              OldZProtContBornF(i,k) = ZERO
              NewZProtContBornB(i,k) = ZERO
              OldZProtContBornB(i,k) = ZERO
 22         continue
 21       continue
 20     continue
      endif

      if (NSolvents .gt. ZERO) then
        do 30 i=1,NSolvents
          nats = NSvnAtoms(i)
          do 31 j=1,nats
            NewZSvnContBorn(i,j) = ZERO
            OldZSvnContBorn(i,j) = ZERO
            NewZSvnContBornF(i,j) = ZERO
            OldZSvnContBornF(i,j) = ZERO
            NewZSvnContBornB(i,j) = ZERO
            OldZSvnContBornB(i,j) = ZERO
 31       continue
 30     continue
      endif

      return
      end
