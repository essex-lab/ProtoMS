c-------------INCLUDED FUNCTION AS AN INLINED BIT OF CODE-------------
c-------------This is because the function call was too expensive!----
c
c This is a rewrite of nbenergy2.inc to only setup charges for a Coulomb energy
c***************************************************************
c       Work out the non-bonded parameters to use for this pair
c***************************************************************
        if (ipar1.lt.0 .or. ipar2.lt.0) then
c         interaction between variable clj parameters
c
          funcsigsclf = .false.
          funcsigsclb = .false.
          if (ipar1.lt.0) then
            funcichg1 = VariableChg(-ipar1,2)
          else
            funcichg1 = charge(ipar1)
          endif
          if (ipar2.lt.0) then
            funcichg2 = VariableChg(-ipar2,2)
          else
            funcichg2 = charge(ipar2)
          endif
c         *** Reference state non bonded parameters ***
          qiqj = funcichg1*funcichg2*FAC4PIEPS0     
c         Calculate forwards and backwards free energies if there is a difference in 
c         lambda
          if (DifferentLambdas) then
            if (ipar1.lt.0) then
               funcichg1b = VariableChg(-ipar1,4)
               funcichg1f = VariableChg(-ipar1,5)
            else
               funcichg1b = funcichg1
               funcichg1f = funcichg1
            endif
            if (ipar2.lt.0) then
               funcichg2b = VariableChg(-ipar2,4)
               funcichg2f = VariableChg(-ipar2,5)
            else
               funcichg2b = funcichg2
               funcichg2f = funcichg2
            endif 
c           ***Forward/Backward state non bonded parameters ***
            qiqjf = funcichg1f*funcichg2f*FAC4PIEPS0
            qiqjb = funcichg1b*funcichg2b*FAC4PIEPS0
          endif
        else
c         interaction between constant clj parameters
          qiqj = ChgIJ(CLJPair(ipar1,ipar2))
          qiqjf = qiqj
          qiqjb = qiqj
c         this is necessary for the q2 sum
          funcichg1 = charge(ipar1)
          funcichg1b = funcichg1 
          funcichg1f = funcichg1 
          funcichg2 = charge(ipar2)
          funcichg2b = funcichg2 
          funcichg2f = funcichg2 
        endif
