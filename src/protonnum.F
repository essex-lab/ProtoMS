      integer function protonNum(ipar)
      include 'dimensions.inc'
      include 'simparams.inc'
      include 'parameters.inc'
c###############################################
c
c     Return the proton number of simulation
c     parameter ipar. If this is a variable 
c     parameter, then the proton number of 
c     lambda=0 is given if lambda<=0.5,
c     and lambda=1 is given if lambda>0.5
c
c     (C) Christopher Woods, October 2004
c
c###############################################

      integer ipar
      
      protonnum = 0
      
      if (ipar.eq.0 .or. ipar.gt.NCLJTypes) then
        return
      else if (ipar.gt.0) then
        protonNum = ProtonNumDB(CLJType(ipar))
      else if (-ipar .le. NVariableCLJ) then
        if (Lambda.le.0.5) then
          protonNum = VariableProtNum(-ipar,1)
        else
          protonNum = VariableProtNum(-ipar,2)
        endif
      endif
      
      return
      end
