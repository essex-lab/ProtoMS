      subroutine readUB(FILE,line,changemode,endfile,nuby)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'constants.inc'
      include 'splitter.inc'
      include 'parameters.inc'
      include 'enums.inc'
c###################################################
c
c     This routine reads in Urey-Bradley parameters -
c     see readBnd and readParFile
c
c     (C) Christopher Woods 
c
c###################################################

      integer FILE
      character*(*) line
      logical changemode,endfile
      integer nuby,iuby
      integer findUreyBradleyIndex,readInt
      double precision readDouble
      character*2 iatm1,iatm2,iatm3
      integer ipar
      double precision ik0,it0
            
      changemode = .false.
      endfile = .false.
      nuby = 0
      
100   continue
        read(FILE,10,err=9700,end=9800) line
10      format(a300)
                
        if (line(1:1).eq.'#') goto 100

c       do not lowercase line as AMBER types are case-sensitive
        
        call split(line)
        
        if (nwords.lt.2) goto 100
        
        call lowercase(words(1))
        
        if (words(1).eq.'mode') then
          changemode = .true.
          call lowercase(words(2))
          return
        else if (words(1).eq.'par') then
c         read in the parameters, format is
c         'par ipar ik0 it0'
          if (nwords.lt.4) then
            call printLine(WARNING,"Not enough parameters for Urey-Bradley parameter line")
            call printLine(WARNING,line)
            goto 100
          endif
          
c         read in the parameters
          ipar = readInt(words(2))
          ik0 = readDouble(words(3))
          it0 = readDouble(words(4))
          
          if (ipar.le.0) then
            call printLine(WARNING,"Cannot have zero or negative Urey-Bradley parameters!")
            call printLine(WARNING,line)
            goto 100
          else if (ipar.gt.MAXUBYPARAM) then
            write(printstring,*) "Parameter ",ipar," is too large (",MAXUBYPARAM,")"
            call printLine(FATAL,printstring)
            write(printstring,*)"Increase MAXUBYPARAM and recompile"
            call closeProgram(-1,printstring)
            goto 100
          endif
          
          if (GotUByParam(ipar)) then
            write(printstring,20) ipar,UByParam(ipar,1),UByParam(ipar,2),ik0,it0
20          format("Overwriting Urey-Bradley parameter ",i4,"(",f8.3,",",f8.3,") with (",
     .                 f8.3,",",f8.3,")")
            call printLine(WARNING,printstring)
            call printLine(WARNING,line)
          endif
          
          UByParam(ipar,1) = ik0
          UByParam(ipar,2) = it0
          GotUByParam(ipar) = .true.
          
          nuby = nuby + 1
        else if (words(1).eq.'atm') then
c         read in the amber atoms, the format it
c         'atm iatm1 iatm2 iatm3 ipar'
          if (nwords.lt.5) then
            call printLine(WARNING,"Not enough atoms in the Urey-Bradley list")
            call printLine(WARNING,line)
            goto 100
          endif
          
c         read in the atoms and parameter number
          iatm1 = words(2)(1:2)
          iatm2 = words(3)(1:2)
          iatm3 = words(4)(1:2)
          ipar = readInt(words(5))
       
          call lowercase(words(2))
          call lowercase(words(3))
          call lowercase(words(4))
          
c         UBs with atoms called 'dm' must be dummies
          if (words(2)(1:2).eq.'dm' .or. words(3)(1:2).eq.'dm' .or. words(4)(1:2).eq.'dm') then
            if (ipar.ne.0) then
              call printLine(WARNING,"Cannot have non-zero parameter for Urey-Bradley that involves dummies!")
              call printLine(WARNING,line)
            endif
            goto 100
          endif
             
          call sortUreyBradley(iatm1,iatm2,iatm3)
          
          iuby = findUreyBradleyIndex(iatm1,iatm2,iatm3)
          
          if (iuby.ne.0) then
c           overwrite the existing uby parameter
            write(printstring,30) iatm1,iatm2,iatm3,UByParamIndex(iuby),ipar
30          format("Overwriting Urey-Bradley parameter for ",a2,"-",a2,"-",a2,", parameter ",
     .             i4," with parameter ",i4)
            call printLine(WARNING,printstring)
            call printLine(WARNING,line)
          else
c           load the new Urey-Bradley parameter
            iuby = NUByParam + 1
            if (NUByParam.gt.MAXUBYPARAM) then
              write(printstring,*) "Exceeded maximum number of Urey-Bradley parameters (",
     .                             MAXUBYPARAM,")"
              call printLine(FATAL,printstring)
              write(printstring,*)"Increase MAXUBYPARAM and recompile"
              call closeProgram(-1,printstring)
              goto 100
            endif
            
            UByName(iuby,1) = iatm1
            UByName(iuby,2) = iatm2
            NUByParam = iuby            
          endif
          UByParamIndex(NUByParam) = ipar
          nuby = nuby + 1                        
        else
          call printLine(WARNING,"Unrecognised Urey-Bradley line!")
          call printLine(WARNING,line)
        endif

      goto 100
9700  continue
      call printLine(WARNING,"Error reading parameter file (UREY-BRADLEY)!")
9800  continue
      endfile = .true.
      return
      
      end
