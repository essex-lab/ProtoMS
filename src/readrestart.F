      subroutine readRestart(restartfile)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'simparams.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'stacks.inc'
      include 'pressure.inc'
      include 'splitter.inc'
      include 'moves.inc'
      include 'parameters.inc'
c######################################################
c
c     This subroutine reads the restart information
c     and sets the coordinates of all of the atoms 
c     appropriately. It then recalculates all of the 
c     zmatricies, and recalculates the total energy
c
c     (C) Christopher Woods   26 November 2002
c
c######################################################

      integer FILE,itmp,i,j,inats,ires,iat,ires2,gcs
      double precision ix,iy,iz,theta
      character*(*) restartfile
      character*300 tmpnam,tmpstr,line
      integer lineLength,lgth,lgth2
      character*4 inam,inam2
      character*4 iresnam,iresnam2
      logical allok, getLogical
      integer getInt
      parameter(FILE=40)

      lgth = lineLength(restartfile)
      if (lgth.eq.0) return
      
      write(printstring,*) "Opening restart file ",restartfile(1:lgth)
      call prettyLine(INFO,printstring)

c     can only open this restart file for reading if it is not already open for writing!
      call closeFile(restartfile)
            
      open(unit=FILE, file=restartfile, status='old', form='formatted',
     .          err=9600)
     
      rewind(FILE)

      allok = .true.
      nwords = 0
c     first line may contain NProteins= or Completed=
c     need slightly more sophisticated text parsing for this
      read(FILE,"(a300)",err=9700,end=9900) line
      call split(line)
      call lowercase(words(1))
      if (words(1)(1:10).eq."completed=") then
         completed = getLogical(words(2),.true.)
         read(FILE,"(11x,i16)",err=9700,end=9900) itmp
         if (.not.completed) then
c     restarting a simulation that did not finish previously
c     adjust starting point of simulation and print details to INFO            
            startstep = itmp
            write(printstring,*) "Restarting simulation that did not complete previously."
            call printLine(INFO,printstring)
            write(printstring,*) "Equilibration will not be performed and simulation step",
     .           " count will begin at ", startstep
            call printLine(INFO,printstring)

c     don't want to resuse the same random seed from the start as with 
c     frequent restarts this could lead to artifacts
c     update the random seed based on the new starting step so its reproducible
            ranseed = ranseed + startstep
            if (ranseed.lt.0) ranseed = -ranseed
            if (mod(ranseed,2).eq.0) ranseed = ranseed - 1
            if (ranseed.lt.100000) ranseed = ranseed+110000
            write(printstring,"(a, i12)") "Random number seed updated based on restart file to ", ranseed
            call prettyLine(INFO,printstring)
            call sgrnd(ranseed)
         endif
      else
c     assume it must have been an NProteins line, grab the integer
c     value and move on to the next bit of logic.
         itmp = getInt(words(2),-1)
         goto 1000
      endif

      
c     read the number of proteins
      read(FILE,10,ERR=9700,END=9900) itmp
 1000 if (itmp.ne.NProteins) then
        write(printstring,*) "The restart file contains the wrong number ",
     .         "of proteins! - ",itmp," vs. ",NProteins
        call printLine(FATAL,printstring)
        allok = .false.
      endif
      
c     now the number of solutes
      read(FILE,10,ERR=9700,END=9900) itmp
      if (itmp.ne.NSolutes) then
        write(printstring,*) "The restart file contains the wrong number ",
     .        "of solutes! - ",itmp," vs. ",NSolutes
        call printLine(FATAL,printstring)
        allok = .false.
      endif
      
c     now the number of solvents
      read(FILE,10,ERR=9700,END=9900) itmp
      if (itmp.ne.NSolvents) then
        write(printstring,*) "The restart file contains the wrong number ",
     .         "of solvents! - ",itmp," vs. ",NSolvents
        call printLine(FATAL,printstring)
        allok = .false.
      endif
      
c     now the number of GCsolutes
      read(FILE,10,ERR=9700,END=9900) itmp
      if (itmp.ne.NGCSolutes) then
        write(printstring,*) "The restart file contains the wrong number ",
     .        "of GCsolutes! - ",itmp," vs. ",NGCSolutes
        call printLine(FATAL,printstring)
        allok = .false.
      endif

      read(FILE,10,ERR=9700,END=9900) itmp
c       PRINT*, "itmp",itmp
       GCSCavity = itmp
10    format(11x,i8)      
     
      if (.not.allok) then
        write(printstring,*) "There were problems reading in the restart file!"
        call closeProgram(-1,printstring)
      endif

c     now read in the box sizes...
      if (boundary.eq.PERIODICBOUNDARIES) then
        read(FILE,72,ERR=9700,END=9900) BoxDimension(1),BoxDimension(2),BoxDimension(3)
        BoxDimension2(1) = 0.5*BoxDimension(1)
        BoxDimension2(2) = 0.5*BoxDimension(2)
        BoxDimension2(3) = 0.5*BoxDimension(3)
        volume = BoxDimension(1)*BoxDimension(2)*BoxDimension(3)
        write(printstring,*) "Solvent box: (",BoxDimension(1),",",BoxDimension(2),",",
     .                        BoxDimension(3),")"
        call printLine(INFO,printstring)
        write(printstring,*) "Volume of box = ",volume
        call prettyLine(INFO,printstring)
      
        read(FILE,72,ERR=9700,END=9900) BoxOrigin(1),BoxOrigin(2),BoxOrigin(3)
        BoxTop(1) = BoxOrigin(1) + BoxDimension(1)
        BoxTop(2) = BoxOrigin(2) + BoxDimension(2)
        BoxTop(3) = BoxOrigin(3) + BoxDimension(3)
72      format(f25.10,1x,f25.10,1x,f25.10)
      endif
      
c     now loop over all proteins, and read in thier coordinates
      do i=1,NProteins
        write(printstring,*) "Reading in restart information for protein ",i
        call printLine(DETAIL,printstring)
        read(FILE,8,ERR=9700,END=9900) tmpnam
8       format(19x,a300)      
        if (tmpnam.ne.ProteinName(i)) then
          write(printstring,*) "Disagreement of protein name in restart file!"
          lgth = lineLength(tmpnam)
          lgth2 = lineLength(ProteinName(i))
          call printLine(WARNING,printstring)
          write(printstring,*) tmpnam(1:lgth)," vs. ",ProteinName(i)(1:lgth2)
          call printLine(WARNING,printstring)
        endif
        
c       read the number of atoms
        read(FILE,12,ERR=9700,END=9900) inats
12      format(9x,i8)
        if (inats.ne.NZAtoms(i)) then
          write(printstring,*) "Disagreement of number of atoms in protein ",i,
     .    " ",inats," vs. ",NZAtoms(i)
          call closeProgram(-1,printstring)
        endif
        
        do j=1,NZAtoms(i)
          ires = Zid2ResID(i,j,1)
          iat = Zid2ResID(i,j,2)
          iresnam = ResNam(i,ires)
          if (iat.gt.0) then
            inam = SCNam(i,ires,iat)
          else if (iat.lt.0) then
            inam = BBNam(i,ires,-iat)
          else
            inam = "????"
          endif           

          read(FILE,11,ERR=9700,END=9900) inam2,iresnam2,ires2,ix,iy,iz
11        format(18x,a4,1x,a3,1x,i8,1x,2(f14.8,1x),f14.8)
          if (inam2.ne.inam .or. iresnam2.ne.iresnam .or.
     .            ires2.ne.ires) then
             write(printstring,*) "Disagreement of protein atom! ",
     .          inam2,"-",iresnam2,"-",ires2," vs. ",inam,"-",iresnam,
     .          "-",ires
             call printLine(WARNING,printstring)
          endif
          
          ZCoords(i,j,1) = ix
          ZCoords(i,j,2) = iy
          ZCoords(i,j,3) = iz
        enddo
      enddo

c     now loop over the solutes and read in their coordinates
      do i=1,NSolutes
        write(printstring,*) "Reading in restart info for solute ",i
        call printLine(DETAIL,printstring)
        
        read(FILE,20,ERR=9700,END=9900) tmpnam
20      format(18x,a300)
        if (tmpnam.ne.SoluteName(i)) then
          lgth = lineLength(tmpnam)
          lgth2 = lineLength(SoluteName(i))
          write(printstring,*) "Disagreement of solute molecule name! ",
     .          tmpnam(1:lgth)," vs. ",SoluteName(i)(1:lgth2)
          call printLine(WARNING,printstring)
        endif
        
        read(FILE,21,ERR=9700,END=9900) inats
21      format(9x,i8)

        if (inats.ne.NSolAtoms(i)) then
          write(printstring,*) "Disagreement of number of atoms in the solute ",
     .        inats," vs. ",NSolAtoms(i)
          call closeProgram(-1,printstring)
        endif

        read(FILE,"(a100)",ERR=9700,END=9900)tmpstr
        if (tmpstr(1:2).eq."Th") then
          read(tmpstr,77) SoluteTheta(i)
77        format(8x,f8.4)

c         read in the dummies and axes for the solute
          read(FILE,23,ERR=9700,END=9900) (SolDummy(i,1,j),j=1,3) 
        else
          SoluteTheta(i) = 1.0d0

c         read in the dummies and axes for the solute
          read(tmpstr,23) (SolDummy(i,1,j),j=1,3)
        endif
        read(FILE,23,ERR=9700,END=9900) (SolDummy(i,2,j),j=1,3)
        read(FILE,23,ERR=9700,END=9900) (SolDummy(i,3,j),j=1,3)
        read(FILE,23,ERR=9700,END=9900) (SolGeomCentre(i,j),j=1,3)
23      format(6x,3(f14.8,1x))        
        
        do j=1,inats
          read(FILE,22,ERR=9700,END=9900) inam2,iresnam2,ires2,ix,iy,iz
22        format(18x,a4,1x,a3,1x,i8,1x,2(f14.8,1x),f14.8)
          if (inam2.ne.SolAtomNam(i,j) .or. ires2.ne.SolResidue(i,j)
     .        .or. iresnam2.ne.SolResNam(i,SolResidue(i,j))) then
             write(printstring,*) "Disagreement of solute atom! ",
     .        inam2,"-",iresnam2,"-",ires2," vs. ",SolAtomNam(i,j),
     .        "-",SolResNam(i,SolResidue(i,j)),"-",SolResidue(i,j)
             call printLine(WARNING,printstring)
          endif
          
          ZSolCoords(i,j,1) = ix
          ZSolCoords(i,j,2) = iy
          ZSolCoords(i,j,3) = iz
        enddo
      enddo

c     now read in the solvent coordinates
      write(printstring,*) "Reading in restart information for the solvent."
      call printLine(DETAIL,printstring)
      
      do i=1,NSolvents
        read(FILE,30,ERR=9700,END=9900) iresnam2
30      format(19x,a4)
        if (iresnam2.ne.SvnName(i)) then
          write(printstring,*) "Disagreement of solvent name ",i,
     .      " ",iresnam2," vs. ",SvnName(i)
          call printLine(WARNING,printstring)
        endif
        
        read(FILE,31,ERR=9700,END=9900) inats
31      format(9x,i8)
        if (inats.ne.NSvnAtoms(i)) then
          write(printstring,*) "Disagreement of number of solvent atoms of ",
     .  "solvent ",i," ",inats," vs. ",NSvnAtoms(i)
          call closeProgram(-1,printstring)
        endif
        
c       read in the solvent axes
        read(FILE,33,ERR=9700,END=9900) (SvnCentre(i,j),j=1,3)
c        read(FILE,33,ERR=9700,END=9900) (SvnMajAxis(i,j),j=1,3)
c        read(FILE,33,ERR=9700,END=9900) (SvnMinAxis(i,j),j=1,3)
c        read(FILE,33,ERR=9700,END=9900) (SvnPerAxis(i,j),j=1,3)
33      format(6x,3(f14.8,1x))
        
        do j=1,inats
          read(FILE,32,ERR=9700,END=9900) inam2,iresnam2,ix,iy,iz
32        format(18x,a4,1x,a3,1x,2(f14.8,1x),f14.8)
          if (inam2.ne.SvnAtomName(i,j)) then
            write(printstring,*) "Disagreement of solvent atom name ",
     .       "in solvent ",i," ",inam2," vs. ",SvnAtomName(i,j)
            call printLine(WARNING,printstring)
          endif
          
          ZSvnCoords(i,j,1) = ix
          ZSvnCoords(i,j,2) = iy
          ZSvnCoords(i,j,3) = iz
        enddo
      enddo
      
c     now read in the GC-solute coordinates
      write(printstring,*) "Reading in restart information for the GCsolutes."
      call printLine(DETAIL,printstring)
      
      do i=1,NGCSolutes
       read(FILE,40,ERR=9700,END=9900) tmpnam
40      format(19x,a18)
c        if (tmpnam.ne.GCSname(i)) then
c          write(printstring,*) "Disagreement of GCsolute name ",i,
c     .      " ",tmpnam," vs. ",GCSname(i)
c          call printLine(WARNING,printstring)
c        endif

        read(FILE,99,ERR=9700,END=9900) theta
99      format(9x,g20.8)
        thetaparam(i) = theta       

      
        read(FILE,41,ERR=9700,END=9900) inats
41      format(9x,i8)
        if (inats.ne.NGCSAtoms(i)) then
          write(printstring,*) "Disagreement of number of GCsolute atoms of ",
     .  "GCsolute ",i," ",inats," vs. ",NGCSAtoms(i)
          call closeProgram(-1,printstring)
        endif

        read(FILE,43,ERR=9700,END=9900) (GCSCentre(i,j),j=1,3)
43      format(6x,3(f14.8,1x))
        
        do j=1,inats
          read(FILE,42,ERR=9700,END=9900) inam2,iresnam2,ix,iy,iz
42        format(18x,a4,1x,a3,1x,2(f14.8,1x),f14.8)
          if (inam2.ne.GCSAtomName(i,j)) then
            write(printstring,*) "Disagreement of GCsolute atom name ",
     .       "in GCsolute ",i," ",inam2," vs. ",GCSAtomName(i,j)
            call printLine(WARNING,printstring)
          endif
          
          ZGCSCoords(i,j,1) = ix
          ZGCSCoords(i,j,2) = iy
          ZGCSCoords(i,j,3) = iz
        enddo
      enddo


c     we have finished with the restart file so it can be safely closed
      close(FILE)      

c     now that all of the coordinates are loaded, we need to rebuild
c     the zmatricies, and recalculate the total energy...
      do i=1,NProteins
        write(printstring,*) "Rebuilding internal information for protein ",i
        call printLine(DETAIL,printstring)
        call syncProteinZmatrix(i)
      enddo

c     need to do the solute and solvent - including the dummy atoms,
c     and the centres/axes!
      do i=1,NSolutes
        if (.not.IgnoreSolute(i)) then        
c         rebuild the solute from these coordinates        
          call getSoluteZmatrix(i)
          call updateSoluteInternals(i)
          call getSoluteCentreRadius(i)
c         synch the forwards and backwards solutes, if there are any
          if (NSolVariable(i).gt.0) call synchFEPSolutes(i)
        endif
      enddo

c     now recalculate the total energy
      GotInitialTotalEnergy = .false.
      call checkGotTotal

      return
      
9600  continue
c     we get here if there was a problem with opening the restart file
      write(printstring,*) "There was a problem opening the restart file! ",
     .                     "Please check to see if the file exists and is readable."
      call closeProgram(-1,printstring)
      return

9700  continue
c     we get here if there was a problem reading a line of the restart file
      write(printstring,*) "There was a problem reading a line of the restart file. ",
     .                     "Is the file corrupted?"
      call closeProgram(-1,printstring)
      return
      
9900  continue
c     we get here if the restart file ended before we were expecting it to!
      write(printstring,*) "The restart file is shorter than was expected. ",
     .                     "Is the file corrupted?"
      call closeProgram(-1,printstring)
      return

      end
