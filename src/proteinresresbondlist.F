      subroutine proteinResResBondList(iprot)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
c#############################################################
c
c     This subroutine scans protein iprot, and picks
c     out all of the residue-residue bonds and stores
c     them in the res-res array...
c
c     The array value is true if the residues are bonded,
c     and false if they are not.
c
c     (C) Christopher Woods   25 January 2002
c
c#############################################################

      integer iprot,nres
      integer i,j,nbnds,nangs,ndihs,nubys
      integer iat1,iat2,iat3,iat4,ires1,ires2,ires3,ires4

      if (iprot.le.0 .or. iprot.gt.MAXPROTEINS) return

c     ResResBond array already initialised to .false. in 
c     initialiseProtein

c     set residues as being bonded to themselves
      nres = NResidues(iprot)
      do i=1,nres
        ResResBond(iprot,i,i) = .true.
      enddo

c     loop over all bonds
      nbnds = NBonds(iprot)
      if (nbnds.gt.0) then
        do i=1,nbnds
          iat1 = Bond(iprot,i,1)
          iat2 = Bond(iprot,i,2)
          ires1 = Zid2ResID(iprot,iat1,1)
          ires2 = Zid2ResID(iprot,iat2,1)
          if (ires1.ne.ires2) then
            ResResBond(iprot,ires1,ires2) = .true.
            ResResBond(iprot,ires2,ires1) = .true.
          endif
        enddo
      endif
      
c     loop over all angles
      nangs = NAngles(iprot)
      if (nangs.gt.0) then
        do i=1,nangs
          iat1 = Angle(iprot,i,1)
          iat2 = Angle(iprot,i,2)
          iat3 = Angle(iprot,i,3)
          ires1 = Zid2ResID(iprot,iat1,1)
          ires2 = Zid2ResID(iprot,iat2,1)
          ires3 = Zid2ResID(iprot,iat3,1)
          if (ires1.ne.ires2) then
            ResResBond(iprot,ires1,ires2) = .true.
            ResResBond(iprot,ires2,ires1) = .true.
          endif
          if (ires1.ne.ires3) then
            ResResBond(iprot,ires1,ires3) = .true.
            ResResBond(iprot,ires3,ires1) = .true.
          endif
          if (ires2.ne.ires3) then
            ResResBond(iprot,ires2,ires3) = .true.
            ResResBond(iprot,ires3,ires2) = .true.
          endif
        enddo
      endif

c     loop over all ubys
      nubys = NUreyBradleys(iprot)
      if (nubys.gt.0) then
        do i=1,nubys
          iat1 = UreyBradley(iprot,i,1)      
          iat2 = UreyBradley(iprot,i,2)
          iat3 = UreyBradley(iprot,i,3)
          ires1 = Zid2ResID(iprot,iat1,1)
          ires2 = Zid2ResID(iprot,iat2,1)
          ires3 = Zid2ResID(iprot,iat3,1)
          if (ires1.ne.ires2) then
            ResResBond(iprot,ires1,ires2) = .true.
            ResResBond(iprot,ires2,ires1) = .true.
          endif
          if (ires1.ne.ires3) then
            ResResBond(iprot,ires1,ires3) = .true.
            ResResBond(iprot,ires3,ires1) = .true.
          endif
          if (ires2.ne.ires3) then
            ResResBond(iprot,ires2,ires3) = .true.
            ResResBond(iprot,ires3,ires2) = .true.
          endif
        enddo
      endif
                          
c     now loop over all dihedrals
      ndihs = NDihedrals(iprot)
      if (ndihs.gt.0) then
        do i=1,ndihs
          iat1 = Dihedral(iprot,i,1)
          iat2 = Dihedral(iprot,i,2)
          iat3 = Dihedral(iprot,i,3)
          iat4 = Dihedral(iprot,i,4)
          ires1 = Zid2ResID(iprot,iat1,1)
          ires2 = Zid2ResID(iprot,iat2,1)
          ires3 = Zid2ResID(iprot,iat3,1)
          ires4 = Zid2ResID(iprot,iat4,1)
          if (ires1.ne.ires2) then
            ResResBond(iprot,ires1,ires2) = .true.
            ResResBond(iprot,ires2,ires1) = .true.
          endif
          if (ires1.ne.ires3) then
            ResResBond(iprot,ires1,ires3) = .true.
            ResResBond(iprot,ires3,ires1) = .true.
          endif
          if (ires1.ne.ires4) then
            ResResBond(iprot,ires1,ires4) = .true.
            ResResBond(iprot,ires4,ires1) = .true.
          endif
          if (ires2.ne.ires3) then
            ResResBond(iprot,ires2,ires3) = .true.
            ResResBond(iprot,ires3,ires2) = .true.
          endif
          if (ires2.ne.ires4) then
            ResResBond(iprot,ires2,ires4) = .true.
            ResResBond(iprot,ires4,ires2) = .true.
          endif
          if (ires3.ne.ires4) then
            ResResBond(iprot,ires3,ires4) = .true.
            ResResBond(iprot,ires4,ires3) = .true.
          endif
        enddo
      endif
         
c     now that all of the inter-residue bonds are found, index all of the
c     bonded residues (index is sorted in numerical order)
      do i=1,nres-1
        do j=i+1,nres
          if (ResResBond(iprot,i,j)) then
            NResResBonds(iprot) = NResResBonds(iprot) + 1
            if (NResResBonds(iprot).gt.MAXRESRESBONDS) then
              write(printstring,*) "Exceeded maximum number of connected residues! (",
     .                        MAXRESRESBONDS,")"
              call printLine(WARNING,printstring)
              call printLine(WARNING,"Increase MAXRESRESBONDS and recompile")
              NResResBonds(iprot) = NResResBonds(iprot) - 1
              return
            endif
            ResResBondIndex(iprot,NResResBonds(iprot),1) = i
            ResResBondIndex(iprot,NResResBonds(iprot),2) = j
          endif
        enddo
      enddo
            
      return
      
      end
