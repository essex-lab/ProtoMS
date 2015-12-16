
      subroutine readParFile(filename)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'splitter.inc'
      include 'proteins.inc'
      include 'sep.inc'
c############################################################
c
c     This routine reads in the new style free format
c     unified ProtoMS parameter files. These files contain
c     the complete description of the forcefield, e.g.
c     parameters and templates!
c
c     (C) Christopher Woods  September 2004
c
c############################################################

      character*(*) filename
      character*300 line,tempstr,pmshome
      integer lineLength,ilgth
c     counters of the number of parameters loaded
      integer i,nclj,iclj,nbnd,ibnd,nang,iang,ndih,idih
      integer nuby,iuby,ntmp,itmp,novr,iovr
c     whether or not we need to change the read mode
      logical hashome,changemode,endfile
      integer FILE
      parameter(FILE=40)
c     JM
      integer isurf,nsurf,iborn,nborn

      call getenv("PROTOMSHOME",pmshome)
      call split(pmshome)
      hashome = nwords.gt.0

      if (hashome.eqv..True.) then
        call split2(filename,PATHSEP)
        do i=1,nwords2
          if (words2(i)(1:lineLength(words2(i))).eq."$PROTOMSHOME") then
            words(i) = pmshome
          else
            words(i) = words2(i)
          endif
        enddo
        nwords = nwords2
        call join2(1,nwords,PATHSEP,tempstr)
        if (filename(1:1).eq.PATHSEP)  then
          write(filename,"(2A)")PATHSEP,tempstr(1:lineLength(tempstr))
        else
          filename = tempstr
        endif
      endif

      ilgth = lineLength(filename)
      if (ilgth.le.0) then
        call closeProgram(-1,"Need to pass filename to readparfile!")
        return
      endif

c     zero the number of parameters read
      nbnd = 0
      nang = 0
      ndih = 0
      nuby = 0
      nclj = 0
      ntmp = 0
      novr = 0
      ibnd = 0
      iang = 0
      idih = 0
      itmp = 0
      iuby = 0
      iclj = 0
      isurf = 0
      nsurf = 0
      iborn = 0
      nborn = 0
      iovr = 0
      
c     open the file      

      write(printstring,"(2a)") "Reading parameter file ",filename(1:ilgth)
      call prettyLine(INFO,printstring)
            
      open(unit=FILE,file=filename,status='old',err=9600)
      rewind(FILE)
      
c     use a 'do while' loop to read all of the lines of the file
100   continue
c       read a line into 'line', jumping to 9700 if there is a problem,
c       or 9800 if we reach the end of the file
        read(FILE,10,err=9700,end=9800) line
10      format(a300)

c       make sure that the line is in lowercase
        call lowercase(line)
        
c       if the line starts with a '#' then skip it
        if (line(1:1).eq.'#') goto 100
        
c       now split the line into words
        call split(line)

c       make sure that the line actually contains something!
        if (nwords.lt.2) goto 100
        
c       now a 'do while' loop while we read over lines of different modes
110     continue        
              
          changemode = .false.
          endfile = .false.

c         see which mode we are in
          if (words(1) .eq. 'mode') then
c           read in different parts based on which mode
            if (words(2) .eq. 'clj') then
c             read in the charge/LJ parameters
              call readCLJ(FILE,line,changemode,endfile,iclj)
c             increment the number of clj parameters
              nclj = nclj + iclj   
            else if (words(2) .eq. 'bond') then
c             read in the bond parameters
              call readBnd(FILE,line,changemode,endfile,ibnd)
              nbnd = nbnd + ibnd
            else if (words(2) .eq. 'angle') then
c             read in the angle parameters
              call readAng(FILE,line,changemode,endfile,iang)
              nang = nang + iang
            else if (words(2) .eq. 'dihedral') then
c             read in the dihedral parameters
              call readDih(FILE,line,changemode,endfile,idih)
              ndih = ndih + idih
            else if (words(2) .eq. 'ureybradley') then
c             read in the urey-bradley terms
              call readUB(FILE,line,changemode,endfile,iuby)
              nuby = nuby + iuby
            else if (words(2) .eq. 'template') then
c             read in templates
              call readTmpl(FILE,line,changemode,endfile,itmp)
              ntmp = ntmp + itmp
c             JM
            else if (words(2) .eq. 'surface') then
c             read in surface terms
              call readSurface(FILE,line,changemode,endfile,isurf)
              nsurf = nsurf + isurf
            else if (words(2) .eq. 'born') then
c             read in generalised born terms
              call readBorn(FILE,line,changemode,endfile,iborn)
              nborn = nborn + iborn
            else if (words(2) .eq. 'override') then
c             read in the overrides
              call readOverride(FILE,line,changemode,endfile,iovr)
              novr = novr + iovr
            else if (words(2) .eq. 'info') then
c             read in forcefield information
              call readInfo(FILE,line,changemode,endfile)            
            else
              write(printstring,"(2a)") "Unknown mode requested in parameter file: ",line
              call printLine(WARNING,line)
            endif                        
          else
            call printLine(WARNING,"Unrecognised line in parameter file:")
            call printLine(WARNING,line)
          endif
          
c         have we just read in a 'mode' line, and thus need to change mode?
          if (changemode) then 
            goto 110
c         have we just reached the end of the file?
          else if (endfile) then 
            goto 9800
c         otherwise we need to read another line
          else
            goto 100
          endif
c         
c       END OF DO-WHILE loop          
         
c     get here when we reach the end of the file
9800  continue
      
c     write some statistics
      write(printstring,20) ntmp,novr,nclj
20    format("Loaded ",I4," templates, ",I4," override templates, ",I4," clj parameters")
      call printLine(INFO,printstring)
      write(printstring,21) nbnd,nang,ndih,nuby
21    format("and ",I4,"/",I4,"/",I4,"/",I4," bond/angle/dihedral/urey-bradley parameters.")
      call printLine(INFO,printstring)
      write(printstring,22) nsurf,nborn
22    format("and ",I4,"/",I4," surface/born parameters.")
      call printLine(INFO,printstring)      
      close(FILE)
      return

c     get here if there was a problem while reading the file
9700  continue
      write(printstring,"(2a)") "There was a problem reading the parameter file ",filename
      call printLine(WARNING,printstring)
      
      return
      
c     get here if there was a problem while opening the file
9600  continue
      write(printstring,"(2a)") "FATAL - There was a problem opening the parameter file ",filename(1:ilgth)
      call closeProgram(-1,printstring)
      
      return
      
      end
