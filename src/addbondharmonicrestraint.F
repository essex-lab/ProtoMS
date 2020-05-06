
      subroutine addBondHarmonicRestraint(itype1,imol1,iatm1,itype2,
     . imol2,iatm2,krestr,bndrestr)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'extraenergies.inc'
c#############################################
c     
c     Add a bond harmonic restraint
c
c     Julien Michel Jan 2007
c
c#############################################
      integer itype1,imol1,iatm1,itype2,imol2,iatm2
      double precision krestr,bndrestr

      nBndHarmonicRestraints = nBndHarmonicRestraints + 1

      BndHarmConst(nBndHarmonicRestraints) = krestr
      BndHarmEq(nBndHarmonicRestraints) = bndrestr
      
      BndHarmAtom1(nBndHarmonicRestraints,1) = itype1
      BndHarmAtom1(nBndHarmonicRestraints,2) = imol1
      BndHarmAtom1(nBndHarmonicRestraints,3) = iatm1

      BndHarmAtom2(nBndHarmonicRestraints,1) = itype2
      BndHarmAtom2(nBndHarmonicRestraints,2) = imol2
      BndHarmAtom2(nBndHarmonicRestraints,3) = iatm2     

      return
      end
