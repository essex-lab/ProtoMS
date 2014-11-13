      subroutine loadSolvent(filename)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'constants.inc'
      include 'splitter.inc'
      include 'simparams.inc'
      include 'solvents.inc'
      include 'flags.inc'
      include 'enums.inc'
c#############################################################
c
c     This routine loads the solvent from a pdb box...
c
c     (C) Christopher Woods, October 2004
c
c#############################################################

      character*(*) filename
      character*300 line
      integer FILE
      integer isvn,iatm,ioldres,nsvn
      integer atnum,iresnum
      character*4 atnam
      character*4 iresnam
      double precision getDouble
      double precision ix,iy,iz
      double precision iorig(3),itop(3),ipos(3),irad,ik
      logical gotwdth,gotcap,ok
      integer lineLength,lgth,linenum
      integer i
      parameter(FILE=40)
      
      lgth = lineLength(filename)
      
      write(printstring,*) "Opening solvent file ",filename(1:lgth)
      call printLine(INFO,printstring)
      
      open(unit=FILE,file=filename,status="old",form='formatted',err=9600)

      rewind(FILE)
      
      isvn = 0
      nsvn = 0
      iatm = 0
      ioldres = 0
      linenum = 0

      gotwdth = .false.
      gotcap = .false.
      ik = ZERO
      irad = ZERO
      do i=1,3
        ipos(i) = ZERO
        iorig(i) = ZERO
        itop(i) = ZERO
      enddo
      
c     read the file in a do-while loop      
100   continue
        linenum = linenum + 1
        read(FILE,10,err=9700,end=9800) line
10      format(a300)

c       lowercase and split the line
        call lowercase(line)
        call split(line)
        
        if (nwords.lt.1) goto 100
        
c       get the solvent box/ball info from the header line
        if (words(1).eq.'header') then
c         the header line should either give the box dimensions, with format;
c         header box xorig yorig zorig xtop ytop ztop
c         or it should give the solvent cap parameters, with format;
c         header cap xpos ypos zpos radius force-constant

          if (words(2).eq.'box') then
            if (gotcap .or. gotwdth) then
              call printLine(WARNING,"Can only have one solvent boundary line per file!")
              call printLine(WARNING,line)
            else if (nwords.lt.5) then
              call printLine(WARNING,"Could not read solvent box info from line;")
              call printLine(WARNING,line)
              call printLine(WARNING,"Format is; header box xorig yorig zorig xtop ytop ztop")
            else if (nwords.lt.8) then
c             assume that the format is 'header box xdim ydim zdim' where we give
c             the dimensions of the box, that is centered on the origin
              iorig(1) = getDouble(words(3),ZERO)
              iorig(2) = getDouble(words(4),ZERO)
              iorig(3) = getDouble(words(5),ZERO)
              
              do i=1,3
                itop(i) = HALF * iorig(i)
                iorig(i) = -itop(i)
              enddo
              
              gotwdth = .true.
            else
              iorig(1) = getDouble(words(3),ZERO)
              iorig(2) = getDouble(words(4),ZERO)
              iorig(3) = getDouble(words(5),ZERO)
              itop(1) = getDouble(words(6),ZERO)
              itop(2) = getDouble(words(7),ZERO)
              itop(3) = getDouble(words(8),ZERO)
              
              gotwdth = .true.
            endif
          else if (words(2).eq.'cap') then
            if (gotcap .or. gotwdth) then
              call printLine(WARNING,"Can only have one solvent boundary line per file!")
              call printLine(WARNING,line)
            else if (nwords.lt.7) then
              call printLine(WARNING,"Could not read solvent cap info from line;")
              call printLine(WARNING,line)
              call printLine(WARNING,"Format is; header cap xpos ypos zpos radius force-constant")
            else
              ipos(1) = getDouble(words(3),ZERO)
              ipos(2) = getDouble(words(4),ZERO)
              ipos(3) = getDouble(words(5),ZERO)
              irad = getDouble(words(6),ZERO)
              ik = getDouble(words(7),ZERO)
              gotcap = .true.
            endif
          endif
        else if (words(1)(1:3).eq.'ter') then
c         we move on to the next solvent molecule        
          iatm = 0
        else
c         try to read in a solvent atom line
          call readPDBLine(line,ok,atnum,atnam,iresnum,iresnam,ix,iy,iz)
          
c         'ok' is true if this line contained a valid pdb atom
          if (.not.ok) then
             if (words(1).eq.'atom' .or. words(1).eq.'hetatm') then
               write(printstring,*) "Could not understand the following PDB line: ",linenum
               call printLine(WARNING,printstring)
               call printLine(WARNING,line)
             endif
            goto 100
          endif

          if (iatm.eq.0 .or. iresnum.ne.ioldres) then
c           this is the first atom of a new solvent...
            isvn = NSolvents + 1
            if (isvn.gt.MAXSOLVENTS) then
              write(printstring,*) "Exceeded maximum number of solvent molecules, ",MAXSOLVENTS
              call printLine(WARNING,printstring)
              call printLine(WARNING,"Increase MAXSOLVENTS and recompile")
              goto 100
            endif
            
            NSolvents = isvn
            nsvn = nsvn + 1
            SvnName(isvn) = iresnam
            SvnRes(isvn) = iresnum
            ioldres = iresnum
            iatm = 0
          endif

          iatm = iatm + 1
          if (iatm.gt.MAXSOLVENTATOMS) then
            write(printstring,*) "Exceeded maximum number of atoms per solvent molecule (",
     .                           MAXSOLVENTATOMS,")"
            call printLine(WARNING,printstring)
            call printLine(WARNING,"Increase MAXSOLVENTATOMS and recompile")
            goto 100
          endif
          
          NSvnAtoms(isvn) = iatm
          SvnAtomName(isvn,iatm) = atnam
          ZSvnCoords(isvn,iatm,1) = ix
          ZSvnCoords(isvn,iatm,2) = iy
          ZSvnCoords(isvn,iatm,3) = iz
        endif
        
        goto 100
c     end of do-while loop

9600  continue
c     get here if we could not open the file
      lgth = lineLength(filename)
      write(printstring,*) "WARNING - unable to load solvent file ",filename(1:lgth)
      call printLine(WARNING,printstring)
      return

9700  continue
c     got here if there was a problem reading a line
      lgth = lineLength(filename)
      write(printstring,*) "WARNING - problem reading line ",linenum," of the solvent file ",
     .     filename(1:lgth)
      call printLine(WARNING,printstring)
      call printLine(WARNING,line)
      call printLine(WARNING,"Skipping the rest of the file!")
      return

c     get here when hit end of the file
9800  continue
c     -------------------------------

      write(printstring,*) "Have now loaded ",NSolvents," solvent molecules (",nsvn,
     .                     " from this file)"
      call printLine(INFO,printstring)

      close(FILE)

c     cannot use boundaries if there are no solvents!      
      if (NSolvents.eq.0) return

c     make sure that only one of gotwdth and gotcap are on (this should be the case)
      if (gotwdth) gotcap = .false.
      if (gotcap) gotwdth = .false.      
            
c     now see what sort of boundaries the user has requested
      if (boundary.eq.NOBOUNDARIES) then 
c       we do not need to do anything!
        gotwdth = .false.
        gotcap = .false.
      else if (boundary.eq.SPHERICALHARMONIC .and. gotwdth) then
c       disagreement between the user request and solvent contents
        call printLine(WARNING,"You have requested the use of a solvent cap but this solvent file")
        call printLine(WARNING,"contains solvent in a box. This is probably not a good idea...")
        gotwdth = .false.
        gotcap = .false.
      else if (gotcap .and. (boundary.eq.CUBICHARMONIC .or. boundary.eq.PERIODICBOUNDARIES)) then
c       again some more disagreement!
        call printLine(WARNING,"You have requested box-type boundaries but this solvent file")
        call printLine(WARNING,"contains information for a solvent cap. This is not good!")
        gotwdth = .false.
        gotcap = .false.
      endif      

      if (gotwdth) then
        if (.not.SolventBox) then
c         we are setting the box dimension for the first time from the solvent file
          do i=1,3
            boxOrigin(i) = iorig(i)
            boxTop(i) = itop(i)
          enddo
        else 
          call printLine(WARNING,"Modifying existing box dimensions by information in the solvent file")
          call printLine(WARNING,"The resulting box may need some NPT equilibration")
c         modify the existing box dimensions - make sure that we use the largest box
          do i=1,3
            if (boxOrigin(i).gt.iorig(i)) boxOrigin(i) = iorig(i)
            if (boxTop(i).lt.itop(i)) boxTop(i) = itop(i)
          enddo
        endif

c       if the boundary is set to agree with the solvent then turn on periodic boundaries                
        if (boundary.eq.AGREESOLVENT) then
          call printLine(INFO,"Periodic boundaries have been turned on in agreement with solvent file")
          boundary = PERIODICBOUNDARIES
        endif

c       turn on the use of a solvent box
        SolventBox = .true.        
                
      else if (gotcap) then
        if (.not.SolventCap) then
c         we are setting the solvent cap from the file for the first time
          capradius = irad
          capconstant = ik           
          do i=1,3
            cap(i) = ipos(i)
          enddo
        else
          call printLine(WARNING,"Cannot set the solvent cap information from multiple sources!")
          call printLine(WARNING,"Only the first source of solvent cap information will be used.")
        endif
        
        if (boundary.eq.AGREESOLVENT) then
          call printLine(INFO,"Spherical cap boundaries have been turned on in agreement with the solvent file")
          boundary = SPHERICALHARMONIC
        endif
        
        SolventCap = .true.
      endif
        
      return
            
      end
