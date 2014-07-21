
      subroutine assignRotamerLibrary(iprot,ires)
      include 'dimensions.inc'
      include 'rotamer.inc'
      include 'proteins.inc'
      include 'printer.inc'
c#########################################################
c
c     This routine assigns the rotamer library dihedrals
c     to residue ires in protein iprot
c
c     (C) Christopher Woods, November 2005
c
c#########################################################

      integer iprot,ires
      integer itmpl,findTemplate,findAtom
      integer nsc,iat1,iat2,iat3,iat4,iatom
      integer i,nassigned
      logical ok
      character*13 inam1,inam2,inam3,inam4,nameString

      if (iprot.le.0 .or. iprot.gt.NProteins) return
      if (ires.le.0 .or. ires.gt.NResidues(iprot)) return

c     first clear the currently assigned rotamer values
      nsc = NSCAtoms(iprot,ires)
      if (nsc.gt.0) then
        do i=1,nsc
          RotamerID(iprot,ires,i) = 0
        enddo
      endif

      nassigned = 0

c     find the template number of this residue name - this is also
c     the id number of the information for this residue in the rotamer library
      itmpl = findTemplate(ResNam(iprot,ires))

c     save the rotamer ID for this residue
      ResRotamerID(iprot,ires) = itmpl

c     if the ID is 0, then something has gone severely wrong!
      if (itmpl.le.0) then
        write(printstring,*) "Cannot find the template for residue ",ResNam(iprot,ires),
     .        ". This is a severe program bug!"
        call printLine(WARNING,printstring)
        return
      endif

c     does this residue have any rotamer states?
      if (NRotamerStates(itmpl) .le. 0 .or. NRotamerDihs(itmpl) .le. 0) then
c       there are no rotamer states - clear the rotamer information for this residue
        ResRotamerID(iprot,ires) = 0
      else
c       there are some rotamer states. We now need to assign the dihedrals
c       Loop over each dihedral in the rotamer library, and try to find it in 
c       this residue
        do i=1,NRotamerDihs(itmpl)
         
c         get the ID numbers of the four atoms of this dihedral
          iat1 = findAtom(iprot,ires,RotamerDihName(itmpl,i,1))
          iat2 = findAtom(iprot,ires,RotamerDihName(itmpl,i,2))
          iat3 = findAtom(iprot,ires,RotamerDihName(itmpl,i,3))
          iat4 = findAtom(iprot,ires,RotamerDihName(itmpl,i,4))

          ok = .true.
           
          if (iat1 .eq. 0) then 
            write(printstring,*) "Could not find rotamer atom ",RotamerDihName(itmpl,i,1),
     .          " in residue ",iprot," ",ResNum(iprot,ires)," ",ResNam(iprot,ires)
            call printLine(WARNING,printstring)
            ok = .false.
          endif

          if (iat2 .eq. 0) then 
            write(printstring,*) "Could not find rotamer atom ",RotamerDihName(itmpl,i,2),
     .          " in residue ",iprot," ",ResNum(iprot,ires)," ",ResNam(iprot,ires)
            call printLine(WARNING,printstring)
            ok = .false.
          endif

          if (iat3 .eq. 0) then 
            write(printstring,*) "Could not find rotamer atom ",RotamerDihName(itmpl,i,3),
     .          " in residue ",iprot," ",ResNum(iprot,ires)," ",ResNam(iprot,ires)
            call printLine(WARNING,printstring)
            ok = .false.
          endif

          if (iat4 .eq. 0) then 
            write(printstring,*) "Could not find rotamer atom ",RotamerDihName(itmpl,i,4),
     .          " in residue ",iprot," ",ResNum(iprot,ires)," ",ResNam(iprot,ires)
            call printLine(WARNING,printstring)
            ok = .false.
          endif

c         make sure that these are four different atoms!
          if (iat1 .eq. iat2 .or. iat1 .eq. iat3 .or. iat1 .eq. iat4 .or.
     .        iat2 .eq. iat3 .or. iat2 .eq. iat4 .or.
     .        iat3 .eq. iat4) then
            call printLine(WARNING,"All of the atoms in the dihedral must be different!")
            ok = .false.
          endif

          if (ok) then
c           we have found all of the atoms in the rotamer dihedral, but do they correspond
c           to a flexible dihedral?

c           sort by higher number (since zmatricies are defined in numerical order...)
c           *do not sort* if the first atom is the O, and the last is 
c           the N, since these are defined in the 'wrong' order...
            if (Zid2ResID(iprot,iat1,2).ne.-4  .and. Zid2ResID(iprot,iat4,2).ne.-1) then
              call sortZDihedral(iat1,iat2,iat3,iat4)
            endif
            
c           now see if this dihedral corresponds to the zmatrix dihedral for the atom...
            iatom = Zid2ResID(iprot,iat1,2)

c           the first atom of the dihedral must be a sidechain atom
            if (iatom .lt. 0) then
              write(printstring,*) "Cannot have a rotamer dihedral on the backbone atom ",
     .          BBNam(iprot,ires, -iatom)
              call printLine(WARNING,printstring)
              ok = .false.

c           the dihedral must correspond to a zmatrix defined dihedral (so that it can 
c           be moved!)
            else if (ZMatSC(iprot,ires,iatom,1) .ne. iat2 .or.
     .               ZMatSC(iprot,ires,iatom,2) .ne. iat3 .or.
     .               ZMatSC(iprot,ires,iatom,3) .ne. iat4) then
              
              inam1 = nameString(iprot,iat1)
              inam2 = nameString(iprot,iat2)
              inam3 = nameString(iprot,iat3)
              inam4 = nameString(iprot,iat4)

              write(printstring,*) "The rotamer dihedral ",RotamerDihName(itmpl,i,1),"-",
     .         RotamerDihName(itmpl,i,2),"-",RotamerDihName(itmpl,i,3),"-",
     .         RotamerDihName(itmpl,i,4)," is not the same as the zmatrix dihedral ",
     .         inam1,"-",inam2,"-",inam3,"-",inam4
              call printLine(WARNING,printstring)
              ok = .false.

c           this is the right dihedral, but it must also be flexible!
            else if ( .not. ZMovSC(iprot,ires,iatom,3) ) then
              call printLine(WARNING,"This dihedral must be flexible to have an assigned rotamer!")
              ok = .false.

c           make sure that this dihedral is not already part of a rotamer
            else if ( RotamerID(iprot,ires,iatom) .ne. 0 ) then
              inam1 = nameString(iprot,iat1)
              inam2 = nameString(iprot,iat2)
              inam3 = nameString(iprot,iat3)
              inam4 = nameString(iprot,iat4)

              write(printstring,*) "Dihedral ",inam1,"-",inam2,"-",inam3,"-",inam4,
     .             " is already part of a rotamer dihedral!"
              call printLine(WARNING,printstring)
              ok = .false.
            
c           this is the right dihedral, and it is flexible! At last! Everything is ok
            else
              RotamerID(iprot,ires,iatom) = i
              nassigned = nassigned + 1
            endif
          endif

          if (.not. ok) then
            write(printstring,*) "Could not assign rotamer ",RotamerDihName(itmpl,i,1),"-",
     .         RotamerDihName(itmpl,i,2),"-",RotamerDihName(itmpl,i,3),"-",
     .         RotamerDihName(itmpl,i,4)
            call printLine(WARNING,printstring)
          endif
        enddo

        if (nassigned .eq. 0) then
          write(printstring,*) "Could not assign any of the rotamer dihedrals to residue ",
     .                          ResNam(iprot,ires)
          call printLine(WARNING,printstring)
          ResRotamerID(iprot,ires) = 0
        endif

      endif

      write(printstring,*) "Assigned ",nassigned," rotamer dihedrals to residue ",
     .       iprot," ",ResNam(iprot,ires)," ",ResNum(iprot,ires)
      call printLine(DETAIL,printstring)

      return 
      end
