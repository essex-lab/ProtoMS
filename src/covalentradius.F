      double precision function covalentRadius(iel)
      implicit none
c##############################################################
c
c     This function will return the covalent radius for element
c     iel (given at the proton(z) number) or 0 if a dummy, 
c     with data taken from 
c
c     http://chemed.chem.purdue.edu/genchem/topicreview/bp/ch7/size.html#cov
c
c     (C) Christopher Woods    1 March 2002
c
c##############################################################

      integer iel
      double precision radii(20)
      
c     data arranged by periodic table z number - only goes up to 
c     Ca - anything above that gets a radius of 1.5A 
c     (noble gases == ZERO, as they do not covalently bond in protoMC)
      data radii / 0.37,                                           0.00,
     .             1.23, 0.89,       0.88, 0.77, 0.70, 0.66, 0.64, 0.00,
     .             1.57, 1.36,       1.25, 1.17, 1.10, 1.04, 0.99, 0.00,
     .             2.03, 1.74 /
      
      if (iel.le.0) then
c       dummies have either a 0, negative of '99' atom type
        covalentRadius = 0.00
        return
      else if (iel.gt.20) then
c       everything greater than 20 has a radius of 1.5 (since the code can't really cope with 
c       'big' atoms anyway.
        covalentRadius = 1.5
        return
      endif
      
      covalentRadius = radii(iel)
      return
      end
