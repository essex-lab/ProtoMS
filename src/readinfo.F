      subroutine readInfo(FILE,line,changemode,endfile)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'splitter.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'flags.inc'
c############################################################
c
c     This routine reads in forcefield information
c     (e.g. what sort of combining rules to use)
c
c     (C) Christopher Woods, October 2004
c
c############################################################

      integer FILE
      character*(*) line
      logical changemode,endfile
      double precision readDouble
      
      changemode = .false.
      endfile = .false.
      
c     read lines in a do-while loop
100   continue
        read(FILE,10,err=9700,end=9800) line
10      format(a300)

        call lowercase(line)
        
c       skip comments
        if (line(1:1).eq.'#') goto 100
        
        call split(line)
        if (nwords.lt.2) goto 100
        
        if (words(1).eq.'mode') then
          changemode = .true.
          return
        else if (words(1).eq.'ljcombine') then
c         reading in the LJ combining rules for sigma parameter
c         Option is either 'geometric' or 'arithmetic'
          if (words(2).eq.'arithmetic') then
            call prettyLine(INFO,"Switching to arithmetic (Lorentz-Berthelot) LJ combining rules")
            lorentzberthelot = .true.
          else if (words(2).eq.'geometric') then
            call prettyLine(INFO,"Switching to geometric (OPLS) LJ combining rules")
            lorentzberthelot = .false.
          else
            call printLine(WARNING,"Could not interpret LJ combining rules from;")
            call printLine(WARNING,line)
          endif
        else if (words(1).eq.'scl14coul') then
c         read in the 1-4 coulomb scale factor
          OneOverScl14C = readDouble(words(2))
          write(printstring,11) OneOverScl14C
11        format("Using a 1-4 coulomb scaling factor of ",f8.4)
          call printLine(INFO,printstring)
        else if (words(1).eq.'scl14lj') then
c         read in the 1-4 LJ scale factor
          OneOverScl14LJ = readDouble(words(2))
          write(printstring,12) OneOverScl14LJ
12        format("Using a 1-4 LJ scaling factor of ",f8.4)
          call printLine(INFO,printstring)
        else
          call printLine(WARNING,"Could not interpret info line;")
          call printLine(WARNING,line)
        endif
        
        goto 100
c     end of do-while loop
9700  continue
      call closeProgram(-1,"There was a problem reading the parameter file!")      
      return
      
9800  continue
c     we have reached the end of the file
      endfile = .true.
      return
                 
      end
