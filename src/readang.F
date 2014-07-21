      subroutine readAng(FILE,line,changemode,endfile,nang)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'splitter.inc'
      include 'parameters.inc'
      include 'enums.inc'
c###################################################
c
c     This routine reads in angle parameters -
c     see readBnd and readParFile
c
c     (C) Christopher Woods 
c
c###################################################

      integer FILE
      character*(*) line
      logical changemode,endfile
      integer nang
      integer findAngleIndex,readInt
      double precision readDouble
      character*2 iatm1,iatm2,iatm3
      integer ipar,iang
      double precision ik0,it0,deg2Rad,rad2deg
    
      changemode = .false.
      endfile = .false.
      nang = 0
      
100   continue
        read(FILE,10,err=9700,end=9800) line
10      format(a300)
        
c       do not lowercase the line as AMBER types are case-sensitive
        
        if (line(1:1).eq.'#') goto 100
        
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
            call printLine(WARNING,"Not enough parameters for angle parameter line")
            call printLine(WARNING,line)
            goto 100
          endif
          
c         read in the parameters
          ipar = readInt(words(2))
          ik0 = readDouble(words(3))
          it0 = readDouble(words(4))
          
          if (ipar.le.0) then
            call printLine(WARNING,"Cannot have zero or negative angle parameters!")
            call printLine(WARNING,line)
            goto 100
          else if (ipar.gt.MAXANGPARAM) then
            write(printstring,*) "Parameter ",ipar," is too large (",MAXANGPARAM,")"
            call printLine(WARNING,printstring)
            call printLine(WARNING,"Increase MAXANGPARAM and recompile")
            goto 100
          endif
          
          if (GotAngParam(ipar)) then
            write(printstring,20) ipar,AngParam(ipar,1),rad2deg(AngParam(ipar,2)),ik0,it0
20          format("Overwriting angle parameter ",i4,"(",f8.3,",",f8.3,") with (",
     .                 f8.3,",",f8.3,")")
            call printLine(WARNING,printstring)
            call printLine(WARNING,line)
          endif
          
          AngParam(ipar,1) = ik0
          AngParam(ipar,2) = deg2Rad(it0)
          GotAngParam(ipar) = .true.
          
          nang = nang + 1
        else if (words(1).eq.'atm') then
c         read in the amber atoms, the format it
c         'atm iatm1 iatm2 iatm3 ipar'
          if (nwords.lt.5) then
            call printLine(WARNING,"Not enough atoms in the angle list")
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
          
c         angles with atoms called 'dm' must be dummies
          if (words(2)(1:2).eq.'dm' .or. words(3)(1:2).eq.'dm' .or. words(4)(1:2).eq.'dm') then
            if (ipar.ne.0) then
              call printLine(WARNING,"Cannot have non-zero angle parameter for angle that involves dummies!")
              call printLine(WARNING,line)
            endif
            goto 100
          endif
             
          call sortAngle(iatm1,iatm2,iatm3)
          
          iang = findAngleIndex(iatm1,iatm2,iatm3)
          
          if (iang.ne.0) then
c           overwrite the existing angle parameter
            write(printstring,30) iatm1,iatm2,iatm3,AngParamIndex(iang),ipar
30          format("Overwriting angle parameter for ",a2,"-",a2,"-",a2,", parameter ",
     .             i4," with parameter ",i4)
            call printLine(WARNING,printstring)
            call printLine(WARNING,line)
          else
c           load the new angle parameter
            iang = NAngParam + 1
            if (iang.gt.MAXANGPARAM) then
              write(printstring,*) "Exceeded maximum number of angle parameters (",
     .                             MAXANGPARAM,")"
              call printLine(WARNING,printstring)
              call printLine(WARNING,"Increase MAXANGPARAM and recompile")
              goto 100
            endif
            
            AngName(iang,1) = iatm1
            AngName(iang,2) = iatm2
            AngName(iang,3) = iatm3
            
            NAngParam = iang
          endif
          
          if (streamHandle(DETAIL).ge.0) then
            write(printstring,40) iatm1,iatm2,iatm2,ipar
40          format("Setting angle parameter between AMBER atoms ",a,"-",a,"-",a," to ",i5)
            call printLine(DETAIL,printstring)
          endif          

          AngParamIndex(iang) = ipar
          
          nang = nang + 1                        
        else
          call printLine(WARNING,"Unrecognised angle line!")
          call printLine(WARNING,line)
        endif

      goto 100
9700  continue
      call printLine(WARNING,"Error reading parameter file (ANGLE)!")
9800  continue
      endfile = .true.
      return
      
      end
