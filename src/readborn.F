      subroutine readBorn(FILE,line,changemode,endfile,nborn)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'splitter.inc'
      include 'gb.inc'
      include 'printer.inc'
c###################################################################
c     This subroutine reads the filehandle 'FILE' looking for born
c     parameters
c     
c     Julien Michel , November 2004
c####################################################################
      integer FILE
      character*(*) line
      logical changemode,endfile
      integer nborn
      integer ipar
      character*2 iatm
      double precision iborn,iscale
      double precision readDouble
      integer readInt

      changemode = .false.
      endfile = .false.
      nborn = 0

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
c         read in a surface parameter - the format of the line is
c         "par ipar iatm irad itens"
          if (nwords.lt.5) then
            call printLine(WARNING,"Not enough parameters for Born line")
            call printLine(WARNING,line)
            goto 100
          endif         
c         read the parameters from the words          
          ipar = readInt(words(2))
          
c         read the AMBER atom type - not that this *is* case-sensitive
          iatm = words(3)(1:2)
          iborn = readDouble(words(4))
          iscale = readDouble(words(5))

          if (ipar.le.0) then
            call printLine(WARNING,"Cannot have zero or negative surface parameters!")
            call printLine(WARNING,line)
            goto 100
          else if (ipar.gt.MAXGBTYPES) then
            write(printstring,*) "Parameter ",ipar," is too large (",MAXGBTYPES,")"
            call printLine(WARNING,printstring)
            call printLine(WARNING,line)
            call printLine(WARNING,"Increase MAXGBTYPES and then recompile.")
            goto 100
          endif
c         print a warning if we are overwriting atom types
          if (GotGB(ipar,1)) then
            write(printstring,20) ipar,GBType(ipar),iatm
 20         format("Overwriting born parameter ",i4,", atm=",A2,", new atm=",A2)
            call printLine(WARNING,printstring)
            call printLine(WARNING,line)
          endif         
c         save the parameters into their arrays
          GBType(ipar) = iatm
          GBIBorn(ipar) = iborn
          GBScale(ipar) = iscale
          nborn = nborn + 1
          GotGB(ipar,1) = .true.

          if (streamHandle(DETAIL).ge.0) then
            write(printstring,40) iatm,ipar,iborn,iscale
 40         format("Setting born parameter for AMBER atom ",A2,
     &        " to ",i3," iborn=",f8.5,"iscale= ",f8.5)
            call printLine(DETAIL,printstring)
          endif 
        else
c         unrecognised line
          call printLine(WARNING,"Unrecognised Born line!")
          call printLine(WARNING,line)            
        endif

        goto 100
c       end of do-while loop
 9700 continue
      call printLine(WARNING,"Error reading parameter file! (Born)")
 9800 continue
c     we have reached the end of the file
      endfile = .true.
      return
      end

