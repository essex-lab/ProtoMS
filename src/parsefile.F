      subroutine parseFile(filename)
      include 'dimensions.inc'
      include 'splitter.inc'
      include 'keyset.inc'
      include 'printer.inc'
c############################################
c
c     This routine parses the ProtoMC 
c     command file and reads each key
c     and its associated option
c
c     (C) Christopher Woods
c
c############################################
      
      character*(*) filename
      character*300 line
      character*128 key
      integer lineLength,lgth

      integer FILE
      parameter(FILE=40)
 
c     initialise number of keys read
      nkeys = 0
      
      lgth = lineLength(filename)
            
      if (lgth.eq.0) then
        call printLine(INFO,"No input command file...")
        return
      else
        write(printstring,*) "Using input command file ",filename(1:lgth)
        call prettyLine(INFO,printstring)
      endif

      open(unit=FILE,file=filename,form='formatted',status='old',err=9600)
    
c     top of do while loop
100   continue
        read(FILE,10,err=9700,end=9800) line
10      format(a300)

c       skip comments
        if (line(1:1).eq.'#') goto 100

c       split the line into words
        call split(line)
        
c       make sure that there are enough words
        if (nwords.lt.2) goto 100        
        
c       lowercase the first word - this is the key
        call lowercase(words(1))
        
c       get the key - this is the first word
        key = words(1)

c       the value is the rest of the line joined back up
        call join(2,nwords,line)      

        lgth = lineLength(line)
        
c       set the key/value pair
        call setKey(key,line)
      
      goto 100

c     end of do while loop

c     --get here if file ends
9800  continue
c     -----------------------

      close(FILE)

      return
      
9600  continue
c     Get here if there was a problem opening the file
      write(printstring,*) "Failed to open input file ",filename(1:lgth)
      call printLine(WARNING,printstring)
      call printLine(WARNING,"Check that the file is readable and that the path/name is correct")
      return

9700  continue
c     Get here if there was a problem reading the file
      write(printstring,*) "There was a problem reading input file ",
     .                     filename(1:lgth)
      call printLine(WARNING,printstring)
      return
      
      end
