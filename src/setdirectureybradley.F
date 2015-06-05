      integer function setDirectUreyBradley(ipar0,ipar1)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'constants.inc'
      include 'simparams.inc'
c####################################################################
c
c     This function directly sets the variable parameter urey-bradley,
c     based on the parameter numbers of the perturbation 
c     (ipar0 is the parameter at lam=0, ipar1 is at lam=1)
c
c     (C) Christopher Woods
c
c####################################################################

      integer ipar0,ipar1,i
      double precision perturbForcefield
      
      setDirectUreyBradley = 0
      
      if (ipar0.lt.0 .or. ipar1.lt.0) then
        write(printstring,*) "Problem in setDirectVariableUreyBradley! ",
     .         "ipar0 = ",ipar0," - ipar1 = ",ipar1
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Urey Bradley parameters *must* be >= 0")
        return
      else if (ipar0.gt.MAXUBYPARAM .or. ipar1.gt.MAXUBYPARAM) then
        write(printstring,*) "Problem in setDirectVariableUreyBradley! ",
     .        "ipar0 = ",ipar0," - ipar1 = ",ipar1," - MAXUBYPARAM = ",MAXUBYPARAM
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Urey Bradley parameters above MAXUBYPARAM are invalid!")
        return
      else if (ipar0.eq.ipar1) then
c       the parameters are the same so this bond does not vary with lambda
        
        if (ipar0.ne.0) then
          if (.not.GotUByParam(ipar0)) then
c           this parameter has not been loaded!
            write(printstring,10) ipar0
10          format("Cannot find Urey-Bradley parameter ",i5,". Check your parameter file to ensure ",
     .             "that this parameter exists. Will now treat this as a dummy!")
            call printLine(WARNING,printstring)
            return
          endif
        endif
        
        setDirectUreyBradley = ipar0
        return
      endif

      if (ipar0.gt.0) then
        if (.not.GotUByParam(ipar0)) then
          write(printstring,10) ipar0
          call printLine(WARNING,printstring)
          return
        endif
      endif
      
      if (ipar1.gt.0) then
        if (.not.GotUByParam(ipar1)) then
          write(printstring,10) ipar1
          call printLine(WARNING,printstring)
          return
        endif
      endif

c     now assign these parameters to the variable bond array
      NVarUreyBradleyParam = NVarUreyBradleyParam + 1
      if (NVarUreyBradleyParam.gt.MAXVARIABLEPARAM) then
        write(printstring,*) "Exceeded number of variable Urey-Bradley parameters (",
     .                         MAXVARIABLEPARAM,")"
        call printLine(FATAL,printstring)
        write(printstring,*)"Increase MAXVARIABLEPARAM and recompile"
        call closeProgram(-1,printstring)
        NVarUreyBradleyParam = NVarUreyBradleyParam - 1
        return
      endif
      
      do i=1,2
        if (ipar0.gt.0) then
          VariableUByParam(NVarUreyBradleyParam,1,i) = UByParam(ipar0,i)
        else
          VariableUByParam(NVarUreyBradleyParam,1,i) = ZERO
        endif
        if (ipar1.gt.0) then
          VariableUByParam(NVarUreyBradleyParam,3,i) = UByParam(ipar1,i)
        else
          VariableUByParam(NVarUreyBradleyParam,3,i) = ZERO
        endif
      enddo

      do i=1,2      
        VariableUByParam(NVarUreyBradleyParam,2,i) = perturbForceField(
     .                               VariableUByParam(NVarUreyBradleyParam,1,i),
     .                               VariableUByParam(NVarUreyBradleyParam,3,i),
     .                                                    Lambda)
      enddo
      
      do i=1,2
        VariableUByParam(NVarUreyBradleyParam,4,i) = perturbForceField(
     .                                VariableUByParam(NVarUreyBradleyParam,1,i),
     .                                VariableUByParam(NVarUreyBradleyParam,3,i),
     .                                                    lambdaB)
        VariableUByParam(NVarUreyBradleyParam,5,i) = perturbForceField(
     .                                VariableUByParam(NVarUreyBradleyParam,1,i),
     .                                VariableUByParam(NVarUreyBradleyParam,3,i),
     .                                                    lambdaF)
      enddo
      
      setDirectUreyBradley = -NVarUreyBradleyParam
      return
      end
