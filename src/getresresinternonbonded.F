      
      subroutine getResResInterNonBonded(iprot,ires1,ires2)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'flags.inc'
      include 'proteins.inc'
c############################################################
c
c     This function gets the interresidue nonbonded atoms 
c     between bonded residues ires1 and ires2
c
c     (C) Christopher Woods, November 2004
c
c############################################################

      integer iprot,ires1,ires2,ibndindex
      integer i,j,istrt1,istrt2,iend1,iend2,n
      integer getResResBondIndex
      integer checkBonded
      integer bondtyp
      character*4 atomName
      
c     get the index number of this bond
      ibndindex = getResResBondIndex(iprot,ires1,ires2)
      if (ibndindex.eq.0) then
        write(printstring,*) "Could not find index of bonded residues ",iprot,ires1,ires2
        call printLine(WARNING,printstring)
        return
      endif      
      
c     the res1 connectivity map is already in TempRes1Map - now we just need
c     to get one for res1
      call createConnectivityMap(iprot,ires2,2)

c     now construct the residue-residue connectivity map
      call createResResConnectivityMap(iprot,ires1,ires2)

c     zero the number of nb pairs between these residues
      NResResNB(iprot,ibndindex) = 0
      NResResNB14(iprot,ibndindex) = 0      
      
      istrt1 = BBZid(iprot,ires1,1)
      iend1 = BBZid(iprot,ires1,4)
      istrt2 = BBZid(iprot,ires2,1)
      iend2 = BBZid(iprot,ires2,4)
            
c     loop over all atoms pairs between residues and see if they are a nbpair, or 14pair.
c     If they are then save them in the ResResNB or ResResNB14 arrays.
      do i=istrt1,iend1
        do j=istrt2,iend2
          bondtyp = checkBonded(iprot,i,j)
          if (bondtyp.gt.0) then
c           these atoms are bonded! - to save space we only record the bonded pairs of
c           atoms between residues (this is different to everything else!). The 
c           function bondedResResEnergy will create the inverse of this array)
            n = NResResNB(iprot,ibndindex) + 1
            if (n.gt.MAXRESRESNB) then
              write(printstring,*) "Exceeded maximum number of bonded-residue bonded pairs. (",
     .                             MAXRESRESNB,")"
              call printLine(FATAL,printstring)
              write(printstring,*)"Increase MAXRESRESNB and recompile"
              call closeProgram(-1,printstring)
            else
              ResResNB(iprot,ibndindex,n,1) = i
              ResResNB(iprot,ibndindex,n,2) = j
              NResResNB(iprot,ibndindex) = n
            endif
          endif
          
          if (bondtyp.eq.3) then
c           these atoms are 1-4 bonded!
            n = NResResNB14(iprot,ibndindex) + 1
            if (n.gt.MAXRESRESNB14) then
              write(printstring,*) "Exceeded maximum number of bonded residue 14-pairs. (",
     .                             MAXRESRESNB14,")"
              call printLine(FATAL,printstring)
              write(printstring,*)"Increase MAXRESRESNB14 and recompile"
              call closeProgram(-1,printstring)
            else
              ResResNB14(iprot,ibndindex,n,1) = i
              ResResNB14(iprot,ibndindex,n,2) = j
              NResResNB14(iprot,ibndindex) = n
            endif
          endif
        enddo
      enddo
      
      if (protodebug) then
        write(printstring,*) "Bonding between ",ResNam(iprot,ires1),"-",ResNum(iprot,ires1)," and ",
     .                        ResNam(iprot,ires2),"-",ResNum(iprot,ires2)
        call printLine(DEBUG,printstring)
        write (printstring,10) (atomName(iprot,i),i=istrt2,iend2)
10      format(5x,50(a4,1x))
        call printLine(DEBUG,printstring)
        
        do i=istrt1,iend1
          write(printstring,11) atomName(iprot,i),(checkBonded(iprot,i,j),j=istrt2,iend2)
11        format(a4,1x,50(1x,i2,2x))
          call printLine(DEBUG,printstring)
        enddo 
        
        write(printstring,12) NResResNB(iprot,ibndindex),NResResNB14(iprot,ibndindex)
12      format("There are ",i4," non-bonded pairs and ",i4," 1-4 pairs.")
        call printLine(DEBUG,printstring)
      endif
              
      return
      end
