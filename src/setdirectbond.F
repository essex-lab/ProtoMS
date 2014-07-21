      integer function setDirectBond(ipar0,ipar1)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'parameters.inc'
      include 'simparams.inc'
      include 'constants.inc'
c####################################################################
c
c     This function directly sets the variable parameter bond,
c     based on the parameter numbers of the perturbation 
c     (ipar0 is the parameter at lam=0, ipar1 is at lam=1)
c
c     (C) Christopher Woods
c
c####################################################################

      integer ipar0,ipar1,i
      double precision perturbForcefield
      
      setDirectBond = 0
      
      if (ipar0.lt.0 .or. ipar1.lt.0) then
        write(printstring,*) "Problem in setDirectBond! ",
     .         "ipar0 = ",ipar0," | ipar1 = ",ipar1
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Bond parameters *must* be >= 0")
        return
      else if (ipar0.gt.MAXBNDPARAM .or. ipar1.gt.MAXBNDPARAM) then
        write(printstring,*) "Problem in setDirectBond! ",
     .        "ipar0 = ",ipar0," | ipar1 = ",ipar1," | MAXBNDPARAM = ",MAXBNDPARAM
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Bond parameters above MAXBNDPARAM are invalid!")
        return
      else if (ipar0.eq.ipar1) then
c       the parameters are the same so this bond does not vary with lambda
        
        if (ipar0.ne.0) then
          if (.not.GotBndParam(ipar0)) then
c           this parameter has not been loaded!
            write(printstring,10) ipar0
10          format("Cannot find bond parameter ",i5,". Check your parameter file to ensure ",
     .             "that this parameter exists. Will now treat this as a dummy!")
            call printLine(WARNING,printstring)
            return
          endif
        endif
        
        setDirectBond = ipar0
        return
      endif

c     make sure that we have loaded the parameters
      if (ipar0.gt.0) then
        if (.not.GotBndParam(ipar0)) then
          write(printstring,10) ipar0
          call printLine(WARNING,printstring)
          return
        endif
      endif
      
      if (ipar1.gt.0) then
        if (.not.GotBndParam(ipar1)) then
          write(printstring,10) ipar1
          call printLine(WARNING,printstring)
          return
        endif
      endif
      
c     now assign these parameters to the variable bond array
      NVarBondParam = NVarBondParam + 1
      if (NVarBondParam.gt.MAXVARIABLEPARAM) then
        write(printstring,*) "Exceeded number of variable bond parameters (",
     .                         MAXVARIABLEPARAM,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Increase MAXVARIABLEPARAM and recompile")
        NVarBondParam = NVarBondParam - 1
        return
      endif
      
      do i=1,2
        if (ipar0.gt.0) then
          VariableBndParam(NVarBondParam,1,i) = BndParam(ipar0,i)
        else
          VariableBndParam(NVarBondParam,1,i) = ZERO
        endif
        if (ipar1.gt.0) then
          VariableBndParam(NVarBondParam,3,i) = BndParam(ipar1,i)
        else
          VariableBndParam(NVarBondParam,3,i) = ZERO
        endif
      enddo

      do i=1,2      
        VariableBndParam(NVarBondParam,2,i) = perturbForceField(
     .                               VariableBndParam(NVarBondParam,1,i),
     .                               VariableBndParam(NVarBondParam,3,i),
     .                                                    Lambda)
      enddo
      
      do i=1,2
        VariableBndParam(NVarBondParam,4,i) = perturbForceField(
     .                                VariableBndParam(NVarBondParam,1,i),
     .                                VariableBndParam(NVarBondParam,3,i),
     .                                                    lambdaB)
      VariableBndParam(NVarBondParam,5,i) = perturbForceField(
     .                                VariableBndParam(NVarBondParam,1,i),
     .                                VariableBndParam(NVarBondParam,3,i),
     .                                                    lambdaF)
      enddo
      
      setDirectBond = -NVarBondParam
      return
      end
