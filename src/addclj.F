      integer function addCLJ(ipar)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'constants.inc'
c##############################################
c
c     This subroutine adds CLJ atom type ipar
c     to the list of CLJ atom types and returns
c     the index of this parameter in the Chg and
c     Sig/Eps arrays.
c
c     (C) Christopher Woods
c
c##############################################

      integer ipar,i

      addCLJ = 0
      
c     if ipar = 0 then it is a dummy atom with no charge or bonding
c     interactions. If it is less than 0, it is a variable parameter
c     and is dealt with elsewhere...
      if (ipar.le.0) return

      if (ipar.gt.MAXCLJ) then
        write(printstring,10) ipar,MAXCLJ
10      format("Requesting a CLJ type(",i4,") above MAXCLJ(",i4,")")
        call printLine(FATAL,printstring)
        write(printstring,*)"Increase MAXCLJ and recompile"
        call closeProgram(-1,printstring)
        return
      endif

c     check to see if this parameter has been loaded...
      if (.not.GotCLJ(ipar,1)) then
        write(printstring,*) "Unknown CLJ parameter ",ipar," requested!"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Make sure that this CLJ parameter is listed in the CLJ parameter section")
        call printLine(WARNING,
     .       "This atom will now be treated as a dummy atom! This is probably not what you want!")
        
c       add this atom as a dummy atom to prevent repeat warnings!
        ChgDB(ipar) = ZERO
        AtmDB(ipar) = "dm"
        SigDB(ipar) = ZERO
        EpsDB(ipar) = ZERO
        ProtonNumDB(ipar) = 0
        GotCLJ(ipar,1) = .true.
      endif
      
c     if this parameter has not already been added, then
      if (.not.GotCLJ(ipar,2)) then
c       if too many have been added then quit.
        NCLJTypes = NCLJTypes + 1
        if (NCLJTypes.gt.MAXCLJTYPES) then
          write(printstring,11) MAXCLJTYPES
11        format("Exceeded maximum number of CLJ types per simulation(",i4,")")
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXCLJTYPES and recompile"
          call closeProgram(-1,printstring)
          return
        endif

        write(printstring,*) "Adding CLJ parameter ",ipar," to the simulation ",ChgDB(ipar),
     .            ", ",AtmDB(ipar),", ",SigDB(ipar),", ",EpsDB(ipar),", ",ProtonNumDB(ipar)
        call printLine(DETAIL,printstring)
        
c       store this type...
        GotCLJ(ipar,2) = .true.
        CLJType(NCLJTypes) = ipar
        addCLJ = NCLJTypes
      else
c       otherwise find the index of the previously loaded parameter
        do i=1,NCLJTypes
          if (ipar.eq.CLJType(i)) then
            addCLJ = i
            return
          endif
        enddo
        
c       we should never get here
        write(printstring,*) "Serious bug in addCLJ - looking for ",ipar
        call closeProgram(-1,printstring)
        addCLJ = 0        

      endif

      return
      
      end
