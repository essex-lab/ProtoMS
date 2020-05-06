      subroutine createConnectivityMap(iprot,ires,map)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'flags.inc'
      include 'proteins.inc'
c##########################################################
c
c     This routine fills the connectivity map 'map'
c     for residue ires in protein iprot
c
c     Atoms are only connected if there is an explicit
c     non-dummy bond between them
c
c     (C) Christopher Woods, October 2004
c
c##########################################################

      integer iprot,ires
c     can't work out how to pass a 2D array!      
      integer map,iat1,iat2
      integer i,j,ibnd,nbnds
      character*4 atomName
      
c     atoms are bonded only if there is an explicit, non-dummy bond between
c     them!
      if (map.eq.1) then
        do i=1,MAXSCATOMS+5
          do j=1,MAXSCATOMS+5
            TempRes1Map(i,j) = .false.
          enddo
        enddo
      else if (map.eq.2) then
        do i=1,MAXSCATOMS+5
          do j=1,MAXSCATOMS+5
            TempRes2Map(i,j) = .false.
          enddo
        enddo
      endif        
        
c     now loop through all of the bonds of this residue and add them to the array
      nbnds = NResBonds(iprot,ires)
      do i=1,nbnds
        ibnd = ResidueBond(iprot,ires,i)
        
c       only work with non-dummy, intra-residue bonds        
        if (.not.BondDummy(iprot,ibnd) .and. Bond(iprot,ibnd,4).eq.ires 
     .         .and. Bond(iprot,ibnd,5).eq.ires) then
c         get the residue based atom numbers for the atoms in this bond
c         - add 5 to allow for the bbatoms, this bbatoms 1-4 are in indicies
c           4 to 1, while sidechain atoms 1+ are in indicies 5+
          iat1 = Zid2ResID(iprot,Bond(iprot,ibnd,1),2) + 5
          iat2 = Zid2ResID(iprot,Bond(iprot,ibnd,2),2) + 5
          
          if (map.eq.1) then
            TempRes1Map(iat1,iat2) = .true.
            TempRes1Map(iat2,iat1) = .true.
          else if (map.eq.2) then
            TempRes2Map(iat2,iat1) = .true.
            TempRes2Map(iat1,iat2) = .true.
          endif
        endif
      enddo

      if (protodebug) then
c       print out the bond map
        write(printstring,*) "Bond map for residue ",ResNam(iprot,ires),"-",ResNum(iprot,ires)
        call printLine(DEBUG,printstring)
        
        write(printstring,10) (atomName(iprot,BBZid(iprot,ires,i)),i=4,1,-1),"????",
     .                        (atomName(iprot,SCZid(iprot,ires,i)),i=1,NSCAtoms(iprot,ires))
10      format(5x,1000(a4,1x))    
        call printLine(DEBUG,printstring)
        
        do i=4,1,-1
          if (map.eq.1) then
            write(printstring,11) atomName(iprot,BBZid(iprot,ires,i)),
     .                     (TempRes1Map(5-i,j),j=1,NSCAtoms(iprot,ires)+5)
          else
            write(printstring,11) atomName(iprot,BBZid(iprot,ires,i)),
     .                     (TempRes2Map(5-i,j),j=1,NSCAtoms(iprot,ires)+5)
          endif
          call printLine(DEBUG,printstring)
        enddo
        
        if (map.eq.1) then
          write(printstring,11) "????",(TempRes1Map(5,j),j=1,NSCAtoms(iprot,ires)+5)
        else
          write(printstring,11) "????",(TempRes2Map(5,j),j=1,NSCAtoms(iprot,ires)+5)
        endif
        call printLine(DEBUG,printstring)
        
11      format(a4,1x,1000(1x,i2,2x))

        do i=1,NSCAtoms(iprot,ires)
          if (map.eq.1) then
            write(printstring,11) atomName(iprot,SCZid(iprot,ires,i)),
     .                     (TempRes1Map(i+5,j),j=1,NSCAtoms(iprot,ires)+5)
          else
            write(printstring,11) atomName(iprot,SCZid(iprot,ires,i)),
     .                     (TempRes2Map(i+5,j),j=1,NSCAtoms(iprot,ires)+5)
          endif
          call printLine(DEBUG,printstring)
        enddo
      endif
              
      return
      end
