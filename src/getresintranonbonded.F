      
      subroutine getResIntraNonBonded(iprot,ires)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'flags.inc'
      include 'proteins.inc'
c############################################################
c
c     This routine gets all of the non-bonded pairs within
c     the residue ires in protein iprot
c
c     (C) Christopher Woods, October 2004
c
c############################################################

      integer iprot,ires
      integer i,j,istrt,iend,n
      integer checkBonded
      integer bondtyp
      character*4 atomName

c     build the connectivity map for this residue in TempRes1Map
      call createConnectivityMap(iprot,ires,1)
            
c     now loop over all atom pairs in the residue
      istrt = BBZid(iprot,ires,1)
      iend = BBZid(iprot,ires,4)
      
c     zero the number of nb pairs in this residue
      NResNB(iprot,ires) = 0
      NResNB14(iprot,ires) = 0
      
      do i=istrt,iend-1
        do j=i+1,iend
c         get the type of bond between atoms i and j
          bondtyp = checkBonded(iprot,i,j)
          if (bondtyp.eq.0) then
c           these atom are not bonded
            n = NResNB(iprot,ires) + 1
            if (n.gt.MAXRESNB) then
              write(printstring,*) "Exceeded maximum number of intra-residue nonbonded pairsxf (",
     .                             MAXRESNB,")"
              call printLine(FATAL,printstring)
              write(printstring,*)"Increase MAXRESNB and recompile"
              call closeProgram(-1,printstring)
            else
              ResNB(iprot,ires,n,1) = i
              ResNB(iprot,ires,n,2) = j
              NResNB(iprot,ires) = n
            endif
          else if (bondtyp.eq.3) then
c           these atoms are 1-4 bonded
            n = NResNB14(iprot,ires) + 1
            if (n.gt.MAXRESNB14) then
              write(printstring,*) "Exceeded maximum number of intra-residue 1-4 pairs (",
     .                             MAXRESNB14,")"
              call printLine(FATAL,printstring)
              write(printstring,*)"Increase MAXRESNB14 and recompile"
              call closeProgram(-1,printstring)
            else
              ResNB14(iprot,ires,n,1) = i
              ResNB14(iprot,ires,n,2) = j
              NResNB14(iprot,ires) = n
            endif
          endif
        enddo
      enddo
      
      if (protodebug) then
        write(printstring,*) "Bonding within ",ResNam(iprot,ires),"-",ResNum(iprot,ires)
        call printLine(DEBUG,printstring)
      
        write(printstring,10) (atomName(iprot,i),i=istrt,iend)
10      format(5x,50(a4,1x))
        call printLine(DEBUG,printstring)
        
        do i=istrt,iend
          write(printstring,11) atomName(iprot,i),
     .                          (checkBonded(iprot,i,j),j=istrt,iend)
11        format(a4,1x,50(1x,i2,2x))
          call printLine(DEBUG,printstring)
        enddo
        write(printstring,12) NResNB(iprot,ires),NResNB14(iprot,ires)
12      format("There are ",i4," non-bonded pairs and ",i4," 1-4 pairs")
        call printLine(DEBUG,printstring)
      endif
        
      return
      end
