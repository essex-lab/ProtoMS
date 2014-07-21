      integer function getDihedral(zpar1,zpar2,zpar3,zpar4)
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'printer.inc'
c####################################################################
c
c     This function sets a variable parameter dihedral
c     with parameter types zat1,zat2,zat3,zat4
c
c     (C) Christopher Woods
c
c#####################################################################

      integer zpar1,zpar2,zpar3,zpar4
      integer ipar0,ipar1,findDihedralParam
      integer setDirectDihedral
      character*2 atA1,atA2,atB1,atB2,atC1,atC2,atD1,atD2
      character*2 amber
      
c#################################################
c#################################################

      getDihedral = 0
      
      if (zpar1.ge.0 .and. zpar2.ge.0 .and. zpar3.ge.0 .and. zpar4.ge.0) then
c       all of the parameters are not variable, so the dihedral is not
c       variable
        atA1 = amber(zpar1)
        atB1 = amber(zpar2)
        atC1 = amber(zpar3)
        atD1 = amber(zpar4)
        getDihedral = findDihedralParam(atA1,atB1,atC1,atD1)
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
      if (zpar4.lt.0) then
        atD1 = VariableAtm(-zpar4,1)
        atD2 = VariableAtm(-zpar4,2)
      else
        atD1 = amber(zpar4)
        atD2 = amber(zpar4)
      endif
      
      ipar0 = findDihedralParam(atA1,atB1,atC1,atD1)
      ipar1 = findDihedralParam(atA2,atB2,atC2,atD2)
c     JM mod
c     If the angle involves one dummy atom, then the parameters are made
c     equal. If you want to change the geometry of the perturbed state 
c     with lambda, then you must set the parameter of the dihedral to to 0 0 
      if (ipar0 .eq. 0) ipar0 = ipar1
      if (ipar1 .eq. 0) ipar1 = ipar0
c      write(printstring,*) 'PARAM FOR ',atA1,'-',atB1,'-',atC1,'-',atD1,
c     & ' to ',atA2,'-',atB2,'-',atC2,'-',atD2,' == ',ipar0,'#',ipar1
c      call printLine(DEBUG,printstring)

c
c     End of JM mod
c
      
      getDihedral = setDirectDihedral(ipar0,ipar1)
      
      return
      
      end
