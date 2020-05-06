      
      subroutine readTmpl(FILE,line,changemode,endfile,ntmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'splitter.inc'
c################################################################
c
c     This routine controls the reading of templates from
c     the parameter file held in filehandle FILE.
c
c     (C) Christopher Woods   September 2004
c
c################################################################

      integer FILE
      character*(*) line
      character*4 resnm
      character*300 solnm
      logical changemode,endfile
      integer ntmpl
      
      changemode = .false.
      endfile = .false.
      ntmpl = 0
      
100   continue
        read(FILE,10,err=9700,end=9800) line
10      format(a300)

        if (line(1:1).eq.'#') goto 100
        
c       lowercase the line
        call lowercase(line)

c       now split into words
        call split(line)
        
        if (nwords.lt.2) goto 100

        if (words(1).eq.'mode') then
          changemode = .true.
          return
        else        
c         the start of yet another fake do-while loop
200       continue
            if (words(1).eq.'residue') then
c             read in a residue template - the name of the residue is the second word (3 letters max)
              resnm = words(2)
              call readResidueTemplate(resnm,FILE,line,changemode,endfile)
              ntmpl = ntmpl + 1
            else if (words(1).eq.'solute') then
c             read in a solute template - the rest of the header line contains the solute name
              call join(2,nwords,solnm)
              call readSoluteTemplate(solnm,FILE,line,changemode,endfile)
              ntmpl = ntmpl + 1
            else if (words(1).eq.'solvent') then
c             read in a solvent template - the name of the solvent is the second word
              resnm = words(2)
              call readSolventTemplate(resnm,FILE,line,changemode,endfile)
              ntmpl = ntmpl + 1
            else if (words(1).eq.'grand') then
             resnm = words(2)
              call readGCsoluteTemplate(resnm,FILE,line,changemode,endfile)
              ntmpl = ntmpl + 1
            else if (words(1).eq.'chain') then
c             read in a backbone template - the name of the backbone is the second word
              solnm = words(2)
              call readChainTemplate(solnm,FILE,line,changemode,endfile)
              ntmpl = ntmpl + 1
            else
              call printLine(WARNING,"Unrecognised template line")
              call printLine(WARNING,line)
              goto 100
            endif    
          
c           if we are changing mode or ending the file, then return
            if (changemode .or. endfile) return
            
c         end of do-while loop
          goto 200
        endif
        
c     end of outer do-while loop
      goto 100
9700  continue
      call printLine(WARNING,"Error reading parameter file (TEMPLATE)")
9800  continue
      endfile = .true.
      return
      
      end
