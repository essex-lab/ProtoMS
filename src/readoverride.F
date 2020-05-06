
      subroutine readOverride(FILE,line,changemode,endfile,novr)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'splitter.inc'
c######################################################################
c
c     This routine reads in the override parameters that are used
c     to re-assign forcefield parameters to already loaded
c     residues. This enables the forcefield to be separated from
c     the structure.
c
c     (C) Christopher Woods, November 2004
c
c######################################################################

      integer FILE
      character*(*) line
      logical changemode,endfile
      integer novr,ipos,interpretPosition
      character*4 resnm
      
      changemode = .false.
      endfile = .false.
      novr = 0
      
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
200       continue
            if (words(1).eq.'residue') then
c             read in a residue override template - format is 'residue name position'
              if (nwords.lt.3) then
c               assume that the position is 'all'
                ipos = ALLPOS
              else 
                ipos = interpretPosition(words(3))              
              endif
              
              resnm = words(2)(1:4)
              
              call readResidueOverride(resnm,ipos,FILE,line,changemode,endfile)
              novr = novr + 1
            else
c             nothing else is currently supported!
              call printLine(WARNING,"Only residue override templates are currently supported!")
              call printLine(WARNING,"Unrecognised residue override template line!")
              call printLine(WARNING,line)
              goto 100
            endif
            
            if (changemode .or. endfile) return
            
          goto 200            
        endif
      goto 100
      
9700  continue
      call printLine(WARNING,"Error reading parameter file (OVERRIDE)")
9800  continue
      endfile = .true.
      
      return
      end
