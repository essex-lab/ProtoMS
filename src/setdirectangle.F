      integer function setDirectAngle(ipar0,ipar1)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'parameters.inc'
      include 'simparams.inc'
      include 'constants.inc'
c####################################################################
c
c     This function directly sets the parameters of a
c     variable parameter angle - see setDirectVariableBond for
c     more explanation
c
c     (C) Christopher Woods   28 September 2002
c
c####################################################################

      integer ipar0,ipar1,i
      double precision perturbForceField
      
      setDirectAngle = 0
      
      if (ipar0.lt.0 .or. ipar1.lt.0) then
        write(printstring,*) "Problem in setDirectAngle! ",
     .         "ipar0 = ",ipar0," - ipar1 = ",ipar1
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Angle parameters *must* be >= 0")
        return
      else if (ipar0.gt.MAXANGPARAM .or. ipar1.gt.MAXANGPARAM) then
        write(printstring,*) "Problem in setDirectAngle! ",
     .        "ipar0 = ",ipar0," | ipar1 = ",ipar1," | MAXANGPARAM = ",MAXANGPARAM
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Angle parameters above MAXANGPARAM are invalid!")
        return
      else if (ipar0.eq.ipar1) then
c       the parameters are the same so this bond does not vary with lambda
        
        if (ipar0.ne.0) then
          if (.not.GotAngParam(ipar0)) then
c           this parameter has not been loaded!
            write(printstring,10) ipar0
10          format("Cannot find angle parameter ",i5,". Check your parameter file to ensure ",
     .             "that this parameter exists. Will now treat this as a dummy!")
            call printLine(WARNING,printstring)
            return
          endif
        endif
        
        setDirectAngle = ipar0
        return
      endif

      if (ipar0.gt.0) then
        if (.not.GotAngParam(ipar0)) then
          write(printstring,10) ipar0
          call printLine(WARNING,printstring)
          return
        endif
      endif
     
      if (ipar1.gt.0) then 
        if (.not.GotAngParam(ipar1)) then
          write(printstring,10) ipar1
          call printLine(WARNING,printstring)
          return
        endif      
      endif
      
c     now add this variable angle
      NVarAngleParam = NVarAngleParam + 1
      if (NVarAngleParam.gt.MAXVARIABLEPARAM) then
        write(printstring,*) "Exceeded number of variable angle parameters (",
     .               MAXVARIABLEPARAM,")"
        call printLine(FATAL,printstring)
        write(printstring,*)"Increase MAXVARIABLEPARAM and recompile"
        call closeProgram(-1,printstring)
        NVarAngleParam = NVarAngleParam - 1
        return
      endif
      
      do i=1,2
        if (ipar0.gt.0) then
          VariableAngParam(NVarAngleParam,1,i) = AngParam(ipar0,i)
        else
          VariableAngParam(NVarAngleParam,1,i) = ZERO
        endif
        if (ipar1.gt.0) then
          VariableAngParam(NVarAngleParam,3,i) = AngParam(ipar1,i)
        else
          VariableAngParam(NVarAngleParam,3,i) = ZERO
        endif
      enddo

      do i=1,2
        VariableAngParam(NVarAngleParam,2,i) = perturbForceField(
     .                                   VariableAngParam(NVarAngleParam,1,i),
     .                                   VariableAngParam(NVarAngleParam,3,i),
     .                                           Lambda)
      enddo

      do i=1,2
          VariableAngParam(NVarAngleParam,4,i) = perturbForceField(
     .                                   VariableAngParam(NVarAngleParam,1,i),
     .                                   VariableAngParam(NVarAngleParam,3,i),
     .                                           LambdaB)
          VariableAngParam(NVarAngleParam,5,i) = perturbForceField(
     .                                   VariableAngParam(NVarAngleParam,1,i),
     .                                   VariableAngParam(NVarAngleParam,3,i),
     .                                           LambdaF)
      enddo
        
      setDirectAngle = -NVarAngleParam
      return
      end
