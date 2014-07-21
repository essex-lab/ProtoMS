
      subroutine readRotamer(filename)
      include 'dimensions.inc'
      include 'rotamer.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'proteins.inc'
      include 'splitter.inc'
c##################################################
c
c     This function reads in the rotamer library
c
c     (C) Christopher Woods, October 2005
c
c##################################################

      character*(*) filename
      character*300 line
      character*4 iresnam
      integer lineLength, ilgth
      integer i,j
      integer nrots
      integer itmpl, findTemplate
      integer iprot,ires,nres
      double precision readDouble, deg2Rad
      integer FILE
      parameter(FILE=42)

c     clear any currently loaded rotamer library
      do i=1,MAXTEMPLATES
        NRotamerStates(i) = 0
        NRotamerDihs(i) = 0
      enddo

      ilgth = lineLength(filename)
      if (ilgth.le.0) then
        call printLine(WARNING,"Need to specify a rotamer file to load!")
        return
      endif

c     open the rotamer library file
      open(UNIT=FILE, file=filename, STATUS="old", err=9600)

c     keep a count of how many rotamer residues have been read
      nrots = 0

c     read in the file in a pseudo do-while loop
100   continue      
        read(FILE,10,err=9700,end=9800) line
10      format(a)

        call split(line)

c       the first word should be the residue name
        if (nwords .le. 0) goto 100
      
c       get the residue name and lowercase it
        iresnam = words(1)(1:4)
        itmpl = 0

        call lowercase(iresnam)     

c       now find the template that corresponds to this residue name
        itmpl = findTemplate(iresnam)

        if (itmpl .eq. 0) then
          write(printstring,*) "Cannot find a residue template for rotamer residue ",
     .                         iresnam
          call printLine(INFO,printstring)
      
c         read away this rotamer (until the next blank line)
200       continue
            read(FILE,10,err=9700,end=9800) line
            call split(line)
            if (nwords .le. 0) goto 100
            goto 200
c         end of reading away the rotamer
        endif

c       ok - we have found a valid rotamer, and its associated template
c       we can now read it in! The first three lines are always present and
c       give nothing, then the number of rotamer states, then the number of
c       rotamer dihedrals for the residue
        nrots = nrots + 1

c       read away the first line - it is not used by anything
        read(FILE,10,err=9700,end=9800) line
        
c       now read in the number of rotamer states
        read(FILE,*,err=9700,end=9800) NRotamerStates(itmpl)

        if (NRotamerStates(itmpl) .gt. MAXROTAMERSTATES) then
          write(printstring,*) "Exceeded maximum number of rotamer states for residue ",
     .         iresnam," (",NRotamerStates(itmpl),") - The maximum number is ",
     .         MAXROTAMERSTATES,". Increase MAXROTAMERSTATES and recompile."
          call closeProgram(-1,printstring)
        endif

c       now read in the number of rotamer dihedrals
        read(FILE,*,err=9700,end=9800) NRotamerDihs(itmpl)

        if (NRotamerDihs(itmpl) .gt. MAXROTAMERDIHS) then
          write(printstring,*) "Exceeded maximum number of rotamer dihedrals for residue ",
     .         iresnam," (",NRotamerDihs(itmpl),") - The maximum number is ",
     .         MAXROTAMERDIHS,". Increase MAXROTAMERDIHS and recompile."
          call closeProgram(-1,printstring)
        endif

c       we can now read in each of the rotamer dihedrals - these are four words, separated 
c       by spaces, e.g. "atm1 atm2 atm3 atm4"
        if (NRotamerDihs(itmpl) .gt. 0) then
          do i=1,NRotamerDihs(itmpl)
            read(FILE,10,err=9700,end=9800) line
            call split(line)
            if (nwords .lt. 4) then
              call join(1,nwords,printstring)
              call printLine(WARNING,"This does not look like a valid rotamer dihedral line!")
              call printLine(WARNING,printstring)
              goto 9700
            endif

c           copy the rotamer dihedral names into storage (lowercase the atom names)
            do j=1,4
              RotamerDihName(itmpl,i,j) = words(j)(1:4)
              call lowercase(RotamerDihName(itmpl,i,j))
            enddo
          enddo
        endif

c       now read in the rotamer states - these are n floating point values, where n is the
c       number of dihedrals in the rotamer
        if (NRotamerStates(itmpl) .gt. 0) then
          do i=1,NRotamerStates(itmpl)
            read(FILE,10,err=9700,end=9800) line
            call split(line)

            if (nwords .lt. NRotamerDihs(itmpl)) then
              call join(1,nwords,printstring)
              call printLine(WARNING,"This does not look like a valid rotamer state line!")
              call printLine(WARNING,printstring)
              goto 9700
            endif

c           copy the rotamer states into storage
            do j=1, NRotamerDihs(itmpl)
              RotamerVal(itmpl,i,j) = readDouble( words(j) )
            enddo
          enddo
        endif

c       we have now finished reading in this rotamer! - print out some info
        write(printstring,*) "Read rotamer library for residue ",iresnam,
     .              " Number of rotamer states = ",NRotamerStates(itmpl),
     .              " Number of rotamer dihedrals = ",NRotamerDihs(itmpl)
        call printLine(INFO,printstring)

c       add a little detail for debugging
        if (NRotamerDihs(itmpl) .gt. 0) then
          call printLine(DETAIL,"Rotamer dihedrals:")
          do i=1,NRotamerDihs(itmpl)
            write(printstring,*) RotamerDihName(itmpl,i,1),"-",RotamerDihName(itmpl,i,2),"-",
     .                           RotamerDihName(itmpl,i,3),"-",RotamerDihName(itmpl,i,4)
            call printLine(DETAIL,printstring)
          enddo 
        else
          call printLine(DETAIL,"No rotamer dihedrals")
        endif    
    
        if (NRotamerStates(itmpl) .gt. 0) then
          call printLine(DETAIL,"Rotamer states:")
          do i=1,NRotamerStates(itmpl)        
            write(printstring,*) (RotamerVal(itmpl,i,j)," ", j=1,NRotamerDihs(itmpl))
            call printLine(DETAIL,printstring)
          enddo
        endif

c       now convert the rotamer states from degrees to radians
        if (NRotamerStates(itmpl) .gt. 0) then
          do i=1,NRotamerStates(itmpl)
            do j=1,NRotamerDihs(itmpl)
              RotamerVal(itmpl,i,j) = deg2Rad(RotamerVal(itmpl,i,j))
            enddo
          enddo
        endif

        goto 100
9800  continue

      close(FILE)

c     get here when we have finished reading the file
      write(printstring,*) "Successfully read in ",nrots," residue rotamer libraries!"
      call printLine(INFO,printstring)

c     now loop over every residue in every protein and assign the read rotamer library
      call printLine(INFO,"Assigning new rotamer library to all loaded protein residues...")
      if (NProteins.gt.0) then
        do iprot=1,NProteins
          nres = NResidues(iprot)
          if (nres .gt. 0) then
            do ires=1,nres
              call assignRotamerLibrary(iprot,ires)
            enddo
          endif
        enddo
      endif

      return

9600  continue
c     get here if there was an error opening the rotamer library file
      write(printstring,*) "Error opening rotamer file ",filename(1:ilgth)
      call printLine(WARNING,printstring)
      close(FILE)

      return

9700  continue
c     get here if there was an error reading the rotamer library file
      ilgth = lineLength(filename)
      write(printstring,*) "Error reading rotamer file ",filename(1:ilgth)
      call printLine(WARNING,printstring)
      close(FILE)

      return
      end
