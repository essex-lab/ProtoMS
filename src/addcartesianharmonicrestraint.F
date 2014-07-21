
      subroutine addCartesianHarmonicRestraint(itype,imol,iatm,
     .        xrestr,yrestr,zrestr,krestr)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'extraenergies.inc'
c#############################################
c     
c     Add a cartesian harmonic restraint ! 
c
c     Julien Michel February 2005
c
c#############################################
      integer itype,imol,iatm
      double precision xrestr,yrestr,zrestr,krestr

      nHarmonicRestraints = nHarmonicRestraints + 1
      HarmPoint(nHarmonicRestraints,1) = xrestr
      HarmPoint(nHarmonicRestraints,2) = yrestr
      HarmPoint(nHarmonicRestraints,3) = zrestr
      HarmConst(nHarmonicRestraints) = krestr
      HarmAtom(nHarmonicRestraints,1) = itype
      HarmAtom(nHarmonicRestraints,2) = imol
      HarmAtom(nHarmonicRestraints,3) = iatm
      
      return
      end
