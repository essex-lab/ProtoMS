      subroutine loadProtein(pfile)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'splitter.inc'
      include 'constants.inc'
      include 'enums.inc'
c##################################################
c
c     this is the subroutine used to load up the 
c     PDB file 'proteinfile' into protein index
c     'iprotein'. There is one index per protein that
c     can be loaded into this code.
c
c     (C) Christopher Woods, October 2004
c
c##################################################

      character*(*) pfile
      integer iprot,ires,iatm,itmpl
      integer atnum,iresnum,oldresnum
      character*4 iresnam,oldresnam
      character*4 atnam
      double precision ix,iy,iz
      integer lgth,lineLength,findResidue,findTemplate,findTemplateAtom
      logical ok
      integer FILE
      integer linenum
      integer lastres,getLastResidue
      integer ipos
      character*300 line

      parameter(FILE=40)

      iprot = NProteins + 1
      
      if (iprot.gt.MAXPROTEINS) then
        write(printstring,*) "Exceeded maximum number of proteins (",MAXPROTEINS,")"
        call printLine(FATAL,printstring)
        write(printstring,*)"Increase MAXPROTEINS and recompile"
        call closeProgram(-1,printstring)
        return
      endif
      
c     first clear space for the protein
      call initialiseProtein(iprot)

      lgth = lineLength(pfile)
      write(printstring,*) "Opening protein file ",pfile(1:lgth)
      call printLine(INFO,printstring)
      open(unit=FILE,file=pfile,status="old",form='formatted',err=9600)
     
      rewind(FILE)

c     get the residue number of the last residue in the chain
      lastres = getLastResidue(FILE,line)
      
      if (lastres.eq.0) then
c       there are no valid residues in this protein!
        call printLine(WARNING,"Could not find a valid protein within this file!")
        goto 9800
      endif     
      
      linenum = 0
      
      oldresnum = 0
      oldresnam = "   "
      ires = 0
      itmpl = 0
      
c     read in the PDB file in a do-while loop
100   continue
        linenum = linenum + 1
        read(FILE,10,err=9700,end=9800) line
10      format(a300)
                
c       lowercase then split the line
        call lowercase(line)
        call split(line)

c       make sure that there are at least two words on the line        
        if (nwords.lt.2) goto 100
        
        if (words(1).eq.'header') then
c         the rest of the header line contains the protein name
          call join(2,nwords,line)
          lgth = lineLength(line)
          ProteinName(iprot) = line(1:lgth)
          write(printstring,*) "Loading protein ",ProteinName(iprot)(1:lgth)
          call printLine(INFO,printstring)
          if (ProteinName(iprot) .eq. 'scoop') then 
            ProteinScoop = .true.
          endif
        else if (words(1)(1:3).eq.'ter') then
c         we cannot load multiple proteins from a protein file!
          write(printstring,*)
     .      "Cannot load multiple proteins from a PDB file - each chain must be in a separate file with ",
     .      "no 'ter' records. Skipping the rest of the PDB file from line ",linenum
          call printLine(WARNING,printstring)
c         the protein has now ended
          goto 9800  
        else if (words(1).eq.'remark') then
          if (words(2)(1:5).eq.'chunk'.or.words(2)(1:9).eq.'pdbparams') then
            call lowercase(words(2))   
c           Increasing the number of keys to store the information found in the protein pdb
            call join(3,nwords,line)      
            call setKey('pdbparams',line)
          endif
        else
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
          
c         see if we have changed residues since the last line
          if (iresnum.eq.oldresnum) then
            if (iresnam.ne.oldresnam) then
              write(printstring,*) "Error in protein PDB on line ",linenum,
     .                         " - two residues have the same residue number but different names!"
              call printLine(WARNING,printstring)
              call printLine(WARNING,line)
              goto 100
            endif
          else
c           we have changed residue number, so search the loaded residues to see if this
c           one has already been loaded
            ires = findResidue(iprot,iresnum)
            if (ires.eq.0) then
c             we have not yet loaded up this residue
              ires = NResidues(iprot) + 1
              if (ires.gt.MAXRESIDUES) then
                write(printstring,*) "Exceeded maximum number of residues per protein (",
     .                    MAXRESIDUES,") on line ",linenum
                call printLine(FATAL,printstring)
                write(printstring,*)"Increase MAXRESIDUES and recompile"
                call closeProgram(-1,printstring)
                goto 100            
              endif
              
c             now find the template for this residue
              itmpl = findTemplate(iresnam)
              if (itmpl.le.0) then
                write(printstring,*) "Unable to find a template for residue ",iresnam," on line ",linenum
                call printLine(WARNING,printstring)
                call printLine(WARNING,"Check to see if this residue is listed in your parameter files")
                goto 100
              endif
              
              NResidues(iprot) = ires
              ResNum(iprot,ires) = iresnum
              ResNam(iprot,ires) = iresnam

c             work out the position of this residue in the chain
c             If this is the first residue, and is also the last residue then
c             this must be a single residue protein (i.e. a single amino acid!)
              if (ires.eq.1 .and. iresnum.eq.lastres) then
                ipos = SINGLE
              else if (ires.eq.1) then
                ipos = FIRST
              else if (iresnum.eq.lastres) then
                ipos = LAST
              else
                ipos = MIDDLE
              endif
c             If dealing with protein scoop then ipos is always MIDDLE (assumes no Nterm/Cterm)
              if (ProteinScoop) ipos = MIDDLE
              
c             save this as the 'old' residue
              oldresnum = iresnum
              oldresnam = iresnam              
            else if (ResNam(iprot,ires).ne.iresnam) then
              write(printstring,*) "Error in protein PDB on line ",linenum,
     .                     " - two residues have the same residue number but different names!"
              call printLine(WARNING,printstring)
              call printLine(WARNING,line)
              goto 100
            endif
          endif  

c         now find out which atom in the template this atom refers to...
          iatm = findTemplateAtom(itmpl,atnam,ipos)
          if (iatm.eq.0) then
            write(printstring,*) "Cannot find atom ",atnam," in template for ",iresnam,
     .                            " from line ",linenum
            call printLine(WARNING,printstring)
            call printLine(WARNING,
     .       "Ensure that this atom is in the residue or chain template for this residue.")
            goto 100
          endif
          
c         if iatm < 0 then this is a BB atom
          if (iatm.lt.0) then
            iatm = -iatm
            if (iatm.gt.4) then
              write(printstring,*) "Serious problem with ProtoMS BB atom assignment from line ",
     .                             "! This is a bug!"
              call printLine(WARNING,printstring)
              goto 100
            endif
c           save the BB atoms data
            BBNam(iprot,ires,iatm) = atnam
            BBNum(iprot,ires,iatm) = atnum
            BBCoord(iprot,ires,iatm,1) = ix
            BBCoord(iprot,ires,iatm,2) = iy
            BBCoord(iprot,ires,iatm,3) = iz
          else
c           this is a side chain atom
            if (iatm.gt.MAXSCATOMS) then
              write(printstring,*) "Exceeded maximum number of sidechain atoms (",
     .                                    MAXSCATOMS,") on line ",linenum
              call printLine(FATAL,printstring)
              write(printstring,*)"Increase MAXSCATOMS and recompile"
              call closeProgram(-1,printstring)
              goto 100
            endif
            SCNam(iprot,ires,iatm) = atnam
            SCNum(iprot,ires,iatm) = atnum
            SCCoord(iprot,ires,iatm,1) = ix
            SCCoord(iprot,ires,iatm,2) = iy
            SCCoord(iprot,ires,iatm,3) = iz
c           make sure that we know how many side chain atoms there are
            if (iatm.gt.NSCAtoms(iprot,ires)) NSCAtoms(iprot,ires) = iatm
          endif
        endif
        
      goto 100
 
9600  continue
      lgth = lineLength(pfile)
      write(printstring,*) "FATAL - unable to load protein PDB file ",pfile(1:lgth)
      call closeProgram(-1,printstring)
      return
9700  continue
      write(printstring,*) "Error while reading protein file line ",linenum,
     .                     ". Skipping the rest of the PDB file!"
      call printLine(WARNING,printstring)
9800  continue

      lgth = lineLength(ProteinName(iprot))
      if (lgth.le.0) then
c       the file did not contain a header line - name the protein after the filename
        lgth = lineLength(pfile)
        if (lgth.gt.120) then
c         only include the end of the file...          
          ProteinName(iprot) = pfile(lgth-119:lgth)
          lgth = 120
        else
          ProteinName(iprot) = pfile
        endif
      endif
      
      write(printstring,20) ProteinName(iprot)(1:lgth),NResidues(iprot)
20    format("Successfully loaded protein ",a,", that contains ",i4," residues")
      call printLine(INFO,printstring)
      
      close(FILE)
      
      NProteins = iprot
                       
      return
      end
