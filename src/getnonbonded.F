      subroutine getNonBonded(iprot)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'enums.inc'
c#############################################################
c
c     This routine controls the calculation of all of the 
c     residue non-bonded pairs in the protein iprot...
c
c     (C) Christopher Woods   13 March 2002
c
c#############################################################

      integer iprot,i,j,nres
      
      if (iprot.le.0 .or. iprot.gt.NProteins) return
      nres = NResidues(iprot)
      if (nres.le.0) return

      write(printstring,*) "Getting the non-bonded list for protein ",iprot
      call printLine(INFO,printstring)      
      do i=1,nres
c       call routine to get intraresidue nb interactions - this sets
c       up the connectivity of ires1 into TempRes1Map
        call getResIntraNonBonded(iprot,i)
        
        do j=i,nres
          if (j.gt.i) then
            if (ResResBond(iprot,i,j)) then
c             these two residues are bonded, so get the interresidue
c             non-bonded pair list
c             - this sets up the connectivity of ires2 into TempRes2Map
c               and the ires1-ires2 connectivity into TempResResMap
              call getResResInterNonBonded(iprot,i,j)
            endif
          endif
        enddo
      enddo
      write(printstring,*) "Finished getting non-bonded list for protein ",iprot
      call printLine(INFO,printstring)
      return
      end
