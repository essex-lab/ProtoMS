      integer function setDirectDihedral(ipar0,ipar1)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'proteins.inc'
c#################################################################
c
c     This function directly sets the parameters of a variable
c     parameter dihedral from the parameter values at lambda = 0
c     and lambda = 1. See setDirectVariableBond for more
c     details
c
c     (C) Christopher Woods
c
c#################################################################

      integer ipar0,ipar1

      setDirectDihedral = 0
      
      if (ipar0.lt.0 .or. ipar1.lt.0) then
        write(printstring,*) "Problem in setDirectDihedral! ",
     .         "ipar0 = ",ipar0," - ipar1 = ",ipar1
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Dihedral parameters *must* be >= 0")
        return
      else if (ipar0.gt.MAXDIHPARAM .or. ipar1.gt.MAXDIHPARAM) then
        write(printstring,*) "Problem in setDirectDihedral! ",
     .        "ipar0 = ",ipar0," | ipar1 = ",ipar1," | MAXDIHPARAM = ",MAXDIHPARAM
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Dihedral parameters above MAXDIHPARAM are invalid!")
        return
      else if (ipar0.eq.ipar1) then
c       the parameters are the same so this dihedral does not vary with lambda
        
        if (ipar0.ne.0) then
          if (.not.GotDihParam(ipar0)) then
c           this parameter has not been loaded!
            write(printstring,10) ipar0
10          format("Cannot find dihedral parameter ",i5,". Check your parameter file to ensure ",
     .             "that this parameter exists. Will now treat this as a dummy!")
            call printLine(WARNING,printstring)
            return
          endif
        endif
        
        setDirectDihedral = ipar0
        return
      endif

      if (ipar0.gt.0) then
        if (.not.GotDihParam(ipar0)) then
          write(printstring,10) ipar0
          call printLine(WARNING,printstring)
          return
        endif
      endif
      
      if (ipar1.gt.0) then
        if (.not.GotDihParam(ipar1)) then
          write(printstring,10) ipar1
          call printLine(WARNING,printstring)
          return
        endif
      endif
          
c     now add this variable dihedral angle
      NVarDihedralParam = NVarDihedralParam + 1
      if (NVarDihedralParam.gt.MAXVARIABLEPARAM) then
        write(printstring,*) "Exceeded number of variable dihedral angle parameters (",
     .                      MAXVARIABLEPARAM,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Increase MAXVARIABLEPARAM and recompile")
        NVarDihedralParam = NVarDihedralParam - 1
        return
      endif
      
c     just save the parameter indexes as the parameters themselves will be calculated
c     on the fly in dihedralEnergy
      VariableDihParam(NVarDihedralParam,1) = ipar0
      VariableDihParam(NVarDihedralParam,2) = ipar1
        
      setDirectDihedral = -NVarDihedralParam
        
      return
      end
