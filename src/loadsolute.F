      subroutine loadSolute(filename)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'constants.inc'
      include 'splitter.inc'
      include 'solutes.inc'
      include 'enums.inc'
c##################################################
c
c     this is the subroutine used to load up the 
c     pdb file filename into the next solute index.
c     There is one index per solute that
c     can be loaded into this code.
c
c     (C) Christopher Woods, October 2004
c     Samuel Genheden implemented reading of several
c                     solutes from a file, 2012
c
c##################################################

      character*(*) filename
      character*300 line,savenam
      integer isol,itmpl
      integer izid
      integer lgth,lineLength,linenum
      integer atnum,iresnum
      integer findMoleculeTemplate,findMolTemplateAtom
      character*4 iresnam
      character*4 atnam
      double precision ix,iy,iz
      logical ok
      integer FILE
      parameter(FILE = 40)
      
      isol = NSolutes + 1
      if (isol.gt.MAXSOLUTES) then
        write(printstring,*) "Exceeded maximum number of solutes (",MAXSOLUTES,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Increase MAXSOLUTES and recompile")
        return
      endif
      
c     first clear space for the protein
      call initialiseSolute(isol)

      lgth = lineLength(filename)
      write(printstring,*) "Opening solute file ",filename(1:lgth)
      call printLine(INFO,printstring)
      open(unit=FILE,file=filename,status="old",form='formatted',err=9600)

c     save this loaded solute
      NSolutes = isol      

      rewind(FILE)
      linenum = 0
      izid = 0
      itmpl = 0
      
c     read the file in a do-while loop
100   continue
        linenum = linenum + 1
        read (FILE,10,err=9700,end=9800) line
10      format(a300)        

c       lowercase the line and split it into words
        call lowercase(line)
        call split(line)
        
c       make sure that there are at least two words on the line
        if (nwords.lt.2) goto 100
        
c       get the solute name from the 'header' line
        if (words(1).eq.'header') then
c         can only have one 'header' line
          if (itmpl.ne.0) then
            call printLine(WARNING,"Can only have one HEADER line in a solute PDB file!")
            call printLine(WARNING,line)
            goto 100
          endif

c         the rest of the header line contains the solute name
          call join(2,nwords,line)
          lgth = lineLength(line)
c         SG
          savenam = line(1:lgth)
          SoluteName(isol) = line(1:lgth)
          write(printstring,*) "Loading solute '",SoluteName(isol)(1:lgth),"'"
          call printLine(INFO,printstring)
          
c         try to find the solute template for this molecule
          itmpl = findMoleculeTemplate(SoluteName(isol)(1:lgth))          
          
          if (itmpl.le.0) then
            write(printstring,*) "Cannot find a template for solute ",SoluteName(isol)(1:lgth)
            call printLine(WARNING,printstring)
            call printLine(WARNING,"Cannot load a solute unless it has a valid template!")
            call prettyLine(WARNING,"This simulation is now broken!")
            return
          endif        
        else if (words(1)(1:3).eq.'ter') then
c         we cannot load multiple solutes from a single file! This could be easily changed though...
c          write(printstring,*)
c     .      "Cannot load multiple solutes from a PDB file - each solute must be in a separate file with ",
c     .      "no 'ter' records. Skipping the rest of the PDB file from line ",linenum
c          call printLine(WARNING,printstring)
c         the solute has now ended
c          goto 9800 
c         SG: when TER record is found, start a new solute molecule
          isol = NSolutes + 1
          if (isol.gt.MAXSOLUTES) then
            write(printstring,*) "Exceeded maximum number of solutes (",MAXSOLUTES,")"
            call printLine(WARNING,printstring)
            call printLine(WARNING,"Increase MAXSOLUTES and recompile")
            return
          endif
          call initialiseSolute(isol)                    
          NSolutes = isol
          izid = 0
          SoluteName(isol) = savenam(1:lineLength(savenam))
        else
          if (itmpl.eq.0) then
            call printLine(WARNING,"Cannot load a solute atom without a valid solute template!")
            call printLine(WARNING,line)
            goto 100
          endif  
        
c         read the pdb line
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

c         try to find this atom in the solute template
          izid = findMolTemplateAtom(itmpl,atnam,iresnam)
          
          if (izid.le.0) then
            write(printstring,*) "Cannot find atom ",atnam,", residue ",iresnam," in the template ",
     .            "for solute ",SoluteName(isol)
            call printLine(WARNING,printstring)
            call printLine(WARNING,"This atom will be ignored!")
            goto 100
          endif
          
c         save the coordinates of the atom           
          ZSolCoords(isol,izid,1) = ix
          ZSolCoords(isol,izid,2) = iy
          ZSolCoords(isol,izid,3) = iz
c          write(printstring,*) ix,iy,iz
c          call printLine(DEBUG,printstring)
c         save the name of the solute atom           
          if (SolAtomNam(isol,izid).eq.atnam) then
            write(printstring,*) "Something strange going on at line ",linenum," of the solute file!"
            call printLine(WARNING,printstring)
            call printLine(WARNING,"An atom-resname seems to appear twice within the PDB file!")
            call printLine(WARNING,"Each atom in the PDB must have a unique atom name / residue name")
            call printLine(WARNING,line)
          endif
          
          SolAtomNam(isol,izid) = atnam
          SolAtomNum(isol,izid) = atnum
                     
c         output the coordinates is 'detail' printing is turned on
          if (streamhandle(DETAIL).ge.0) then
            write(printstring,11) SolAtomNum(isol,izid),SolAtomNam(isol,izid),iresnam,
     .            ZSolCoords(isol,izid,1),ZSolCoords(isol,izid,2),
     .            ZSolCoords(isol,izid,3)
11          format("ATOM ",i5,1x,a,1x,a,1x,3(1x,f8.3))
            call printLine(DETAIL,printstring) 
          endif
        endif
        goto 100
c     end of do-while loop
        
c     --- get here at the end of the file
9800  continue
c     -----------------------------------

c     make sure that the solute has a name
      if (SoluteName(isol).eq.' ') then
        lgth = lineLength(filename)
        write(printstring,*) "You haven't named the solute in file ",filename(1:lgth)
        call printLine(WARNING,printstring)
        call printLine(WARNING,"You need to name the solute so that ProtoMC can find its template.")
        call printLine(WARNING,"You name a solute by adding the line;")
        call printLine(WARNING,"header solutename")
        call printLine(WARNING,"to the pdb file (where 'solutename' is the name of the solute).")
        call printLine(WARNING,"ProtoMC will now continue, naming the solute after the filename")
        SoluteName(isol) = filename(1:lgth)
      endif
      
c     now that the solute has been read, close the file
      close(FILE)

      return
     
9600  continue
c     get here if we could not open the file
      lgth = lineLength(filename)
      write(printstring,*) "WARNING - unable to load solute file ",filename(1:lgth)
      call printLine(WARNING,printstring)
      return

9700  continue
c     got here if there was a problem reading a line
      lgth = lineLength(filename)
      write(printstring,*) "WARNING - problem reading line ",linenum," of the solute file ",
     .     filename(1:lgth)
      call printLine(WARNING,printstring)
      call printLine(WARNING,line)
      call printLine(WARNING,"Skipping the rest of the file!")
      return
      
      end
