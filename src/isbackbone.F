      logical function isBackBone(iprot,izatm)
      include 'dimensions.inc'
      include 'proteins.inc'
c####################################################
c
c     This is a simple function that returns whether
c     an atom in iprot protein, with zmatrix id
c     izatm, is a backbone atom..
c
c     (C) Christopher Woods   25 January 2002
c
c####################################################

      integer iprot,izatm

      isBackBone = .false.
      
      if (Zid2ResID(iprot,izatm,2).lt.0) isBackBone = .true.
      
      return
      
      end
