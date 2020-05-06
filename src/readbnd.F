      subroutine readBnd(FILE,line,changemode,endfile,nbnd)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'splitter.inc'
      include 'parameters.inc'
      include 'enums.inc'
c##############################################################
c
c     This routine reads in bond parameters from the 
c     FILE filehandle
c
c     (C) Christopher Woods,  September 2004
c
c##############################################################

      integer FILE
      character*(*) line
      logical changemode,endfile
      integer nbnd
      double precision ir0,ik0
      integer ipar,ibnd
      character*2 iatm1,iatm2
      double precision readDouble
      integer readInt,findBondIndex
            
      changemode = .false.
      endfile = .false.
      nbnd = 0
      
100   continue
        read(FILE,10,err=9700,end=9800) line
10      format(a300)
        
c       do not lowercase the line as amber atom types are case-sensitive
        
c       skip comments
        if (line(1:1).eq.'#') goto 100
        
c       split the line into words
        call split(line)
        
c       make sure that the line contains something
        if (nwords.lt.2) goto 100
        
c       lowercase the first words
        call lowercase(words(1))        
        
        if (words(1).eq.'mode') then
          changemode = .true.
c         also need to lowercase the second word
          call lowercase(words(2))          
          return
        else if (words(1).eq.'par') then
c         read in parameter values - format is
c         "par ipar ik0 ir0"
          if (nwords.lt.4) then
            call printLine(WARNING,"Not enough parameters for a bond parameter line")
            call printLine(WARNING,line)
            goto 100
          endif
      
c         read in the parameters
          ipar = readInt(words(2))
          ik0 = readDouble(words(3))
          ir0 = readDouble(words(4))
          
          if (ipar.le.0) then
            call printLine(WARNING,"Cannot have zero or negative bond parameters!")
            call printLine(WARNING,line)
            goto 100
          else if (ipar.gt.MAXBNDPARAM) then
            write(printstring,*) "Parameter ",ipar," is too large (",MAXBNDPARAM,")"
            call printLine(FATAL,printstring)
            write(printstring,*)"Increase MAXBNDPARAM and recompile"
            call closeProgram(-1,printstring)
            goto 100
          endif
          
          if (GotBndParam(ipar)) then
            write(printstring,20) ipar,BndParam(ipar,1),BndParam(ipar,2),ik0,ir0
20          format("Overwriting bond parameter ",i4,"(",f8.3,",",f8.3,") with (",
     .                 f8.3,",",f8.3,")")
            call printLine(WARNING,printstring)
            call printLine(WARNING,line)
          endif
          
          BndParam(ipar,1) = ik0
          BndParam(ipar,2) = ir0
          GotBndParam(ipar) = .true.
          
          nbnd = nbnd + 1
        else if (words(1).eq.'atm') then
c         read in the amber atoms - format is
c         "atm  atm1  atm2 ipar"
          if (nwords.lt.4) then
            call printLine(WARNING,"Not enought atoms in the bond list")
            call printLine(WARNING,line)
            goto 100
          endif
          
c         read in the atoms
          iatm1 = words(2)(1:2)
          iatm2 = words(3)(1:2)
          ipar = readInt(words(4))

c         if either atoms are 'DM' (dummies) then their bond parameter *must* be zero
          call lowercase(words(2))
          call lowercase(words(3))

          if (words(2)(1:2).eq.'dm' .or. words(3)(1:2).eq.'dm') then
            if (ipar.ne.0) then
              call printLine(WARNING,"Cannot have non-zero bond parameter for bond that involves dummies!")
              call printLine(WARNING,line)
            endif
            goto 100
          endif
          
          call sortBond(iatm1,iatm2)

          ibnd = findBondIndex(iatm1,iatm2)
                              
          if (ibnd.ne.0) then
c           overwrite the existing bond parameter
            write(printstring,30) iatm1,iatm2,BndParamIndex(ibnd),ipar
30          format("Overwriting bond parameter for ",a2,"-",a2,", parameter ",
     .             i4," with parameter ",i4)
            call printLine(WARNING,printstring)
            call printLine(WARNING,line)
          else
c           create a new bond parameter
            ibnd = NBndParam + 1
            if (ibnd.gt.MAXBNDPARAM) then
              write(printstring,*) "Exceeded maximum number of bond parameters! (",
     .                              MAXBNDPARAM,")"
              call printLine(FATAL,printstring)
              write(printstring,*)"Increase MAXBNDPARAM and recompile"
              call closeProgram(-1,printstring)
              goto 100
            endif
            
            BndName(ibnd,1) = iatm1
            BndName(ibnd,2) = iatm2
            NBndParam = ibnd
          endif
          
          if (streamHandle(DETAIL).ge.0) then
            write(printstring,40) iatm1,iatm2,ipar
40          format("Setting bond parameter between AMBER atoms ",a,"-",a," to ",i5)
            call printLine(DETAIL,printstring)
          endif          

          BndParamIndex(ibnd) = ipar
          nbnd = nbnd + 1          
        else
          call printLine(WARNING,"Unrecognised bond line!")
          call printLine(WARNING,line)
        endif
        
        goto 100
c     end of do-while loop
9700  continue
      call printLine(WARNING,"Error reading parameter file (BOND)!")
9800  continue
      endfile = .true.
      return
      
      end
