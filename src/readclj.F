      subroutine readCLJ(FILE,line,changemode,endfile,nclj)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'constants.inc'
      include 'splitter.inc'
      include 'parameters.inc'
      include 'enums.inc'
c##########################################################
c
c     This subroutine reads the filehandle 'FILE' looking
c     for clj parameters. 
c
c     (C) Christopher Woods,  September 2004
c
c##########################################################

      integer FILE
      character*(*) line
      logical changemode,endfile
      integer nclj
      double precision ichg,isig,ieps
      character*2 iatm
      double precision readDouble
      integer readInt
      integer ipar,iprot
      
      changemode = .false.
      endfile = .false.
      nclj = 0
      
c     loop reading lines from the file      
100   continue
        read(FILE,10,err=9700,end=9800) line
10      format(a300)

c       do not lowercase the line as the AMBER atom type are case-sensitive
        
c       skip comments
        if (line(1:1).eq.'#') goto 100
        
c       split the line
        call split(line)
c       make sure that the line contains something
        if (nwords.lt.2) goto 100

c       lowercase the first two words
        call lowercase(words(1))     
        call lowercase(words(2))   
                
        if (words(1).eq.'mode') then
c         we are changing mode
          changemode = .true.
          return
        else if (words(1).eq.'par') then
c         read in a clj parameter - the format of the line is
c         "par ipar iatm iprot ichg  isig  ieps"
          if (nwords.lt.7) then
            call printLine(WARNING,"Not enough parameters for CLJ line")
            call printLine(WARNING,line)
            goto 100
          endif
c         read the parameters from the words          
          ipar = readInt(words(2))
          
c         read the AMBER atom type - not that this *is* case-sensitive
          iatm = words(3)(1:2)
          iprot = readInt(words(4))
          ichg = readDouble(words(5))
          isig = readDouble(words(6))
          ieps = readDouble(words(7))

          if (ipar.le.0) then
            call printLine(WARNING,"Cannot have zero or negative CLJ parameters!")
            call printLine(WARNING,line)
            goto 100
          else if (ipar.gt.MAXCLJ) then
            write(printstring,*) "Parameter ",ipar," is too large (",MAXCLJ,")"
            call printLine(FATAL,printstring)
            write(printstring,*)"Increase MAXCLJ and recompile"
            call closeProgram(-1,printstring)
            goto 100
          endif

c         print a warning if we are overwriting atom types
          if (GotCLJ(ipar,1)) then
            write(printstring,20) ipar,AtmDB(ipar),iatm
20          format("Overwriting clj parameter ",i4,", atm=",A2,", new atm=",A2)
            call printLine(WARNING,printstring)
            call printLine(WARNING,line)
          endif
                         
c         save the parameters into their arrays
          ChgDB(ipar) = ichg
          AtmDB(ipar) = iatm
          SigDB(ipar) = isig
          EpsDB(ipar) = ieps
          ProtonNumDB(ipar) = iprot
          GotCLJ(ipar,1) = .true.

          nclj = nclj + 1
        else
c         unrecognised line
          call printLine(WARNING,"Unrecognised CLJ line!")
          call printLine(WARNING,line)
        endif
        
        goto 100
c     end of do-while loop
9700  continue
      call printLine(WARNING,"Error reading parameter file! (CLJ)")
9800  continue
c     we have reached the end of the file
      endfile = .true.
      return
      end
