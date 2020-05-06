      subroutine pushBackPDAcont
c####################################
c     Save the pairwise contribution 
c     to the effective born radii of 
c     each atom
c  
c     Julien Michel November 2004
c####################################
      include 'dimensions.inc'
      include 'constants.inc'
      include 'parameters.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'solvents.inc'      
      include 'gb.inc'
      integer i,j,k,nats,nres,resstart,resend

      if (NSolutes .gt. ZERO) then
        do 10 i=1,NSolutes
          nats = NSolAtoms(i)
          do 20 j=1,nats
            StackZSolContBorn(i,j) = ZSolContBorn(i,j)
            StackZSolContBornF(i,j) = ZSolContBornF(i,j)
            StackZSolContBornB(i,j) = ZSolContBornB(i,j)
 20       continue
 10     continue
      endif

      if (NProteins .gt. ZERO) then
        do 30 i=1,NProteins
          nres = NResidues(i)
          do 40 j=1,nres
            resstart = BBZid(i,j,1)
            resend = BBZid(i,j,4)
            do 50 k=resstart,resend
              StackZProtContBorn(i,k) = ZProtContBorn(i,k)
              StackZProtContBornF(i,k) = ZProtContBornF(i,k)
              StackZProtContBornB(i,k) = ZProtContBornB(i,k)
 50         continue
 40       continue
 30     continue
      endif

      if (NSolvents .gt. ZERO) then
        do 60 i=1,NSolvents
          nats = NSvnAtoms(i)
          do 70 j=1,nats
            StackZSvnContBorn(i,j) = ZSvnContBorn(i,j)
            StackZSvnContBornF(i,j) = ZSvnContBornF(i,j)
            StackZSvnContBornB(i,j) = ZSvnContBornB(i,j)
 70       continue
 60     continue
      endif

      return
      end
