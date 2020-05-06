      subroutine TotalBorn(savenrg)
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
c##############################################################
c This subroutine computes the born radii of the entire system
c     
c     Julien Michel November 2004
c##############################################################
      integer i,j,k,l,nres,nresk,savenrg
      logical hasSoluteMoved,hasResidueMoved,hasSolventMoved
c      integer nats,resstart,resend

c
c  Perform the accumulation of the pairwise contributions
c
      if (NProteins.gt.0) then      
        do 100 i=1,NProteins
          nres = NResidues(i)
          do 101 j=1,nres
c         Do pair pair intra protein    
c            write (printstring,*) 'DOING PROT ',i,' RES INTRA ',j
c           call printLine(DEBUG,printstring)    
            if ( hasResidueMoved(i,j) ) call ResidueBornIntraPDA(savenrg,i,j)
 101      continue
c         Now the contributions between other residues in the same protein
          if (nres .gt. 1) then
            do 102 j=1,nres-1
             do 103 k=j+1,nres
               if ( hasResidueMoved(i,j) .or. hasResidueMoved(i,k) )
     .           call ResidueResidueBornInterPDA(savenrg,i,j,i,k)
c                 call ResidueResidueBornInterPDA(savenrg,i,j,i,k)
c                 write (printstring,*) 'REDOING PROT ',i,' RES ',j,' RES ',k
c                 call printLine(DEBUG,printstring)        
c               endif
 103         continue
 102        continue
          endif
 100    continue
      endif

      if (NProteins .gt. 1) then
c       Do the contributions between a residue and all the residues in every 
c       other protein
        do 104 i=1,NProteins-1
          nres = NResidues(i)
          do 105 j=1,nres
            do 106 k=i+1,NProteins
              nresk = NResidues(k)
              do 107 l=1,nresk
c                write (printstring,*) 'DOING PROT ',i,' RES ',j,
c     &            ' PROT ',k,' RES ',l
c                call printLine(DEBUG,printstring)      
                if ( hasResidueMoved(i,j) .or. hasResidueMoved(k,l) )
     .            call ResidueResidueBornInterPDA(savenrg,i,j,k,l)
 107          continue
 106        continue
 105      continue
 104    continue
      endif

      if (NSolutes.gt.0) then
        do 200 i=1,NSolutes
c         Do pair-pair intrasolute
          if (IgnoreSolute(i)) goto 200
          if ( hasSoluteMoved(i) ) call SoluteBornIntraPDA(savenrg,i)
 200    continue

        if (NSolutes.gt.1) then
          do 201 i=1,NSolutes-1
            if (IgnoreSolute(i)) goto 201
            do 202 j=i+1,NSolutes
c         Do pair-pair inter solutes
c              write (printstring,*) 'DOING SOLUTE SOLUTE ',i,'-',j
c              call printLine(DEBUG,printstring)
              if (IgnoreSolute(j)) goto 202
              if ( hasSoluteMoved(i) .or. hasSoluteMoved(j) ) 
     .          call SoluteSoluteBornInterPDA(savenrg,i,j)
 202        continue
 201      continue
        endif
c
c Now do the solutes proteins contributions
c
        if (NProteins.gt.0) then
          do 210 i=1,NProteins
            nres = NResidues(i)
            do 220 j=1,nres
              do 230 k=1,NSolutes
                if (IgnoreSolute(k)) goto 230
c                write (printstring,*) 'DOING PROT ',i,' RES ',j,
c     &            ' SOLUTE ',k
c                call printLine(DEBUG,printstring)
                if ( hasResidueMoved(i,j) .or. hasSoluteMoved(k) )
     .            call residueSoluteBornPDA(savenrg,i,j,k)
 230          continue
 220        continue
 210      continue
        endif
      endif

c     now do the same for the solvent
      if (NSolvents.gt.0) then
        do 300 i=1,NSolvents
c          write (printstring,*) 'DOING SOLVENT INTRA ',i
c          call printLine(DEBUG,printstring)
          if ( hasSolventMoved(i) ) call SolventBornIntraPDA(savenrg,i)
 300    continue
        
        if (NSolvents .gt. 1) then
          do 301 i=1,NSolvents-1
            do 302 j=i+1,NSolvents
c              write (printstring,*) 'DOING SOLVENT SOLVENT ',i,'-',j
c              call printLine(DEBUG,printstring)
              if ( hasSolventMoved(i) .or. hasSolventMoved(j) )
     .         call SolventSolventBornInterPDA(savenrg,i,j)              
 302        continue
 301      continue
        endif
c       Now the solvent proteins 
c        
        if (NProteins.gt.0) then
          do 310 i=1,NProteins
            nres = NResidues(i)
            do 311 j=1,nres
              do 312 k=1,NSolvents
c                write (printstring,*) 'DOING PROT ',i,' RES ',j,
c     &            ' SOLVENT ',k
c                call printLine(DEBUG,printstring)
                if ( hasResidueMoved(i,j) .or. hasSolventMoved(k) )
     .            call ResidueSolventBornPDA(savenrg,i,j,k)                
 312          continue
 311        continue
 310      continue
        endif
c Do the solute solvent contributions
        if (NSolutes .gt. 0) then
          do 320 i=1,NSolutes
            if (IgnoreSolute(i)) goto 320
            do 330 j=1,NSolvents
c                write (printstring,*) 'DOING SOLUTE ',i,' SOLVENT ',j
c                call printLine(DEBUG,printstring)
              if ( hasSoluteMoved(i) .or. hasSolventMoved(j) )
     .          call SoluteSolventBornPDA(savenrg,i,j)               
 330        continue
 320      continue
        endif
      endif

      return
      end        
