      integer function getUreyBradley(zpar1,zpar2,zpar3)
      include 'dimensions.inc'
      include 'parameters.inc'
c######################################################
c
c     This function adds a variable parameter 
c     urey-bradley term to the simulation between
c     atoms with simulation parameters zpar1 zpar2 and zpar3
c
c     (C) Christopher Woods
c
c######################################################

      integer zpar1,zpar2,zpar3,ipar0,ipar1
      character*2 atA1,atA2,atB1,atB2,atC1,atC2
      character*2 amber
      integer setDirectUreyBradley
      integer findUreyBradleyParam
      
c#################################################
c#################################################
      
      getUreyBradley = 0
      
      if (zpar1.ge.0 .and. zpar1.ge.0 .and. zpar3.ge.0) then
c       all of the parameters are fixed with lambda, so
c       the UBy must also be fixed with lambda
        atA1 = amber(zpar1)
        atB1 = amber(zpar2)
        atC1 = amber(zpar3)
        getUreyBradley = findUreyBradleyParam(atA1,atB1,atC1)
        return
      endif
      
c     get the names of the combinations of atoms
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
      
c     get the bond parameter at lambda = 0
      ipar0 = findUreyBradleyParam(atA1,atB1,atC1)
c     now the parameter at lambda = 1
      ipar1 = findUreyBradleyParam(atA2,atB2,atC2)
c     JM mod
      if (ipar0 .eq. 0) ipar0 = ipar1
      if (ipar1 .eq. 0) ipar1 = ipar0

      getUreyBradley = setDirectUreyBradley(ipar0,ipar1)

      return
      end
