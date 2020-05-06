      subroutine createResResConnectivityMap(iprot,ires1,ires2)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'flags.inc'
      include 'proteins.inc'
c######################################################################
c
c     This routine creates the connectivity map between atoms
c     of residues ires1 and ires2 in protein iprot
c
c     (C) Christopher Woods, October 2004
c
c######################################################################
      
      integer nbnds,i,j,ibnd
      integer iat1,iat2,iprot,ires1,ires2
      character*4 atomName
      
c     atoms are bonded only if there is an explicit, non-dummy bond between
c     them

c     initialise the map      
      do i=1,MAXSCATOMS+5
        do j=1,MAXSCATOMS+5
          TempResResMap(i,j) = .false.
        enddo
      enddo
            
c     scan the bonds belonging to residue 1;
c     atoms in residue 1 are only bonded to those in residue
c     2 if there is an explicit, non-dummy bond linking them
      nbnds = NResBonds(iprot,ires1)
      do i=1,nbnds
        ibnd = ResidueBond(iprot,ires1,i)
        
c       only work with non-dummy, intra-residue bonds        
        if (.not.BondDummy(iprot,ibnd)) then
          if (Bond(iprot,ibnd,4).eq.ires1 .and. 
     .        Bond(iprot,ibnd,5).eq.ires2) then
c           get the residue based atom numbers for the atoms in this bond
c           - add 5 to allow for the bbatoms, this bbatoms 1-4 are in indicies
c             4 to 1, while sidechain atoms 1+ are in indicies 5+
            iat1 = Zid2ResID(iprot,Bond(iprot,ibnd,1),2) + 5
            iat2 = Zid2ResID(iprot,Bond(iprot,ibnd,2),2) + 5
            TempResResMap(iat1,iat2) = .true.
          else if (Bond(iprot,ibnd,4).eq.ires2 .and.
     .             Bond(iprot,ibnd,5).eq.ires1) then
            iat2 = Zid2ResID(iprot,Bond(iprot,ibnd,1),2) + 5
            iat1 = Zid2ResID(iprot,Bond(iprot,ibnd,2),2) + 5
            TempResResMap(iat1,iat2) = .true.
          endif
        endif
      enddo

      if (protodebug) then
        write(printstring,*) "Bondmap between ",ResNam(iprot,ires1),"-",ResNum(iprot,ires1),
     .                       " and ",ResNam(iprot,ires2),"-",ResNum(iprot,ires2)
        call printLine(DEBUG,printstring)
        
        write(printstring,10) (atomName(iprot,BBZid(iprot,ires2,i)),i=4,1,-1),"????",
     .                        (atomName(iprot,SCZid(iprot,ires2,i)),i=1,NSCAtoms(iprot,ires2))
10      format(5x,1000(a4,1x))
        call printLine(DEBUG,printstring)
        
        do i=4,1,-1
          write(printstring,11) atomName(iprot,BBZid(iprot,ires1,i)),
     .                          (TempResResMap(5-i,j),j=1,NSCAtoms(iprot,ires2)+5)
          call printLine(DEBUG,printstring)
        enddo

11      format(a4,1x,1000(1x,i2,2x))        
                
        write(printstring,11) "????",(TempResResMap(5,j),j=1,NSCAtoms(iprot,ires2)+5)
        call printLine(DEBUG,printstring)

        do i=1,NSCAtoms(iprot,ires1)
          write(printstring,11) atomNAme(iprot,SCZid(iprot,ires1,i)),
     .                          (TempResResMap(i+5,j),j=1,NSCAtoms(iprot,ires2)+5)
          call printLine(DEBUG,printstring)
        enddo
      endif
      
      return
      end
