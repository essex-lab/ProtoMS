      integer function getAngle(zpar1,zpar2,zpar3)
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'printer.inc'
c####################################################################
c
c     This function sets a variable parameter angle
c     with simulation parameter types zat1,zat2,zat3
c
c     (C) Christopher Woods
c
c#####################################################################

      integer zpar1,zpar2,zpar3,ipar0,ipar1,findAngleParam
      character*2 atA1,atA2,atB1,atB2,atC1,atC2
      character*2 amber
      integer setDirectAngle
      
c#################################################
c#################################################

      getAngle = 0
      
      if (zpar1.ge.0 .and. zpar2.ge.0 .and. zpar3.ge.0) then
c       none of the parameters vary with lambda, so the
c       angle won't vary with lambda
        atA1 = amber(zpar1)
        atB1 = amber(zpar2)
        atC1 = amber(zpar3)
        getAngle = findAngleParam(atA1,atB1,atC1)
        return
      endif
      
      if (zpar1.lt.0) then
        atA1 = VariableAtm(-zpar1,1)
        atA2 = VariableAtm(-zpar1,2)
      else
        atA1 = amber(zpar1)
        atA2 = amber(zpar1)
      endif
      if (zpar2.lt.0) then
        atB1 = VariableAtm(-zpar2,1)
        atB2 = VariableAtm(-zpar2,2)
      else
        atB1 = amber(zpar2)
        atB2 = amber(zpar2)
      endif
      if (zpar3.lt.0) then
        atC1 = VariableAtm(-zpar3,1)
        atC2 = VariableAtm(-zpar3,2)
      else
        atC1 = amber(zpar3)
        atC2 = amber(zpar3)
      endif
      
c     get the parameters for the angle at lam=0 and lam=1
      ipar0 = findAngleParam(atA1,atB1,atC1)
      ipar1 = findAngleParam(atA2,atB2,atC2)

c
c     JM mod
c     If the angle involves one dummy atom, then the parameters are made
c     equal

      if (ipar0 .eq. 0) ipar0 = ipar1
      if (ipar1 .eq. 0) ipar1 = ipar0
c      write(printstring,*) 'PARAM FOR ',atA1,'-',atB1,'-',atC1,' to ',
c     & atA2,'-',atB2,'-',atC2,' == ',ipar0,'#',ipar1
c      call printLine(DEBUG,printstring)
c
c     End of JM mod
c

      getAngle = setDirectAngle(ipar0,ipar1)
      
      return
      
      end
