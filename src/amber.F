      character*2 function amber(ipar)
      include 'dimensions.inc'
      include 'simparams.inc'
      include 'parameters.inc'
c###############################################
c
c     Return the amber atom name of simulation
c     parameter ipar. If this is a variable 
c     parameter, then the amber name of 
c     lambda=0 is given if lambda<=0.5,
c     and lambda=1 is given if lambda>0.5
c
c     (C) Christopher Woods, October 2004
c
c###############################################

      integer ipar
      
      amber = "dm"
      
      if (ipar.eq.0 .or. ipar.gt.NCLJTypes) then
        return
      else if (ipar.gt.0) then
        amber = AtmDB(CLJType(ipar))
      else if (-ipar .le. NVariableCLJ) then
        if (Lambda.le.0.5) then
          amber = VariableAtm(-ipar,1)
        else
          amber = VariableAtm(-ipar,2)
        endif
      endif
      
      return
      end
