      
      integer function getBond(zpar1,zpar2)
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'printer.inc'
c######################################################
c
c     This function adds bond between atoms with
c     simulation parameters zpar1 and zpar2 to
c     the simulation, and returns the index of 
c     the added bond parameter
c
c     (C) Christopher Woods, October 2004
c
c######################################################

      integer zpar1,zpar2,ipar0,ipar1
      character*2 atA0,atA1,atB0,atB1
      character*2 amber
      integer setDirectBond
      integer findBondParam
      
c#################################################
c#################################################
      
      getBond = 0
      
      if (zpar1.ge.0 .and. zpar2.ge.0) then
c       neither of the parameters vary with lambda, so
c       it is easy to get the amber types and the bond
c       cannot be variable        
        atA1 = amber(zpar1)
        atB1 = amber(zpar2)
        getBond = findBondParam(atA1,atB1)
        return
      endif

c     one of the atom parameters is variable, so the bond parameter may be
c     variable...
            
c     get the names of the combinations of atoms
      if (zpar1.lt.0) then
        atA0 = VariableAtm(-zpar1,1)
        atA1 = VariableAtm(-zpar1,2)
      else
        atA0 = amber(zpar1)
        atA1 = amber(zpar1)
      endif
      
      if (zpar2.lt.0) then
        atB0 = VariableAtm(-zpar2,1)
        atB1 = VariableAtm(-zpar2,2)
      else
        atB0 = amber(zpar2)
        atB1 = amber(zpar2)
      endif
      
c     get the bond parameter at lambda = 0
      ipar0 = findBondParam(atA0,atB0)
c     now the parameter at lambda = 1
      ipar1 = findBondParam(atA1,atB1)


c     JM mod
c     If the bond involves one dummy atom, then the parameters are made
c     equal. If you want to change that dof with lambda, 
c     then you must set the parameter of the bond to 0 0
c     otherwise, be ready for some large free energy changes ! 

      if (ipar0 .eq. 0) ipar0 = ipar1
      if (ipar1 .eq. 0) ipar1 = ipar0

c      write(printstring,*) 'PARAM FOR ',atA0,'-',atB0,' to ',
c     &  atA1,'-',atB1,' == ',ipar0,'#',ipar1
c      call printLine(DEBUG,printstring)
c     End of JM mod
c
      getBond = setDirectBond(ipar0,ipar1)

      return
      end
