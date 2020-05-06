      subroutine updateBornRadii
c##############################################
c Set Old Effective Born = New Effective Born for 
c the atoms whose 
c GB energy has been recomputed
c
c     Julien Michel November 2004
c##############################################
      include 'dimensions.inc'
      include 'constants.inc'
      include 'parameters.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'solvents.inc'      
      include 'gb.inc'
      integer i,j,k,nats,nres,resstart,resend

      if (NSolutes.gt.0) then
        do 10 i=1,NSolutes
          nats = NSolAtoms(i)
          do 20 j=1,nats
            if (SolAtomdoBorn(i,j)) then 
              OldZSolEffBorn(i,j) = ZSolEffBorn(i,j)
              OldZSolEffBornF(i,j) = ZSolEffBornF(i,j)
              OldZSolEffBornB(i,j) = ZSolEffBornB(i,j)
            endif
 20       continue
 10     continue
      endif

      if (NProteins .gt. 0) then
        do 30 i=1,NProteins
          nres = NResidues(i)
          do 40 j=1,nres
            resstart = BBZid(i,j,1)
            resend = BBZid(i,j,4)
            do 50 k=resstart,resend
              if (ProtAtomdoBorn(i,k)) then 
                OldZProtEffBorn(i,k) = ZProtEffBorn(i,k)
                OldZProtEffBornF(i,k) = ZProtEffBornF(i,k)
                OldZProtEffBornB(i,k) = ZProtEffBornB(i,k)
              endif
 50         continue
 40       continue
 30     continue
      endif

      if (NSolvents .gt. ZERO) then
        do 60 i=1,NSolvents
          nats = NSvnAtoms(i)
          do 70 j=1,nats
            if (SvnAtomdoBorn(i,j)) then
              OldZSvnEffBorn(i,j) = ZSvnEffBorn(i,j)
              OldZSvnEffBornF(i,j) = ZSvnEffBornF(i,j)
              OldZSvnEffBornB(i,j) = ZSvnEffBornB(i,j)
            endif
 70       continue
 60     continue
      endif

      return
      end
