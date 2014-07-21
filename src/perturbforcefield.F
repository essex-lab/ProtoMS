      double precision function perturbForceField(val0,val1,lam)
      include 'dimensions.inc'
      include 'constants.inc'
c###############################################################
c
c     This function uses simple linear scaling to morph
c     the forcefield from val0 at lam=0, to val1 at lam=1.
c     The value at lam=vallam is returned.
c
c     (C) Christopher Woods    28 January 2002
c
c###############################################################

      double precision val0,val1,lam

      if (lam.eq.ZERO) then
        perturbForceField = val0
      else if (lam.eq.ONE) then
        perturbForceField = val1
      else
        perturbForceField = val0 + (lam*(val1-val0))
      endif
      
      return
      
      end
