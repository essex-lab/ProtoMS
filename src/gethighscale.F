      double precision function getHighScale(lam)
      implicit none
c############################################################
c
c     Returns the scale factor for the high solute at 
c     lambda=lam. Currently just does linear scaling, but
c     could be modified to do a whole range of other scaling
c
c     The sum of getLowScale and getHighScale for each
c     value of lambda should be equal to one!
c
c     (C) Christopher Woods
c
c############################################################

      double precision lam
      
      getHighScale = lam
      return
      end
