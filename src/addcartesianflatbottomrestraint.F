
      subroutine addCartesianFlatBottomRestraint(itype,imol,iatm,
     .        xrestr,yrestr,zrestr,krestr,wrestr)
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
      double precision xrestr,yrestr,zrestr,krestr,wrestr

      nFlatBottomRestraints = nFlatBottomRestraints + 1
      FBotPoint(nFlatBottomRestraints,1) = xrestr
      FBotPoint(nFlatBottomRestraints,2) = yrestr
      FBotPoint(nFlatBottomRestraints,3) = zrestr
      FBotConst(nFlatBottomRestraints) = krestr
      FBotWidth(nFlatBottomRestraints) = wrestr
      FBotAtom(nFlatBottomRestraints,1) = itype
      FBotAtom(nFlatBottomRestraints,2) = imol
      FBotAtom(nFlatBottomRestraints,3) = iatm
      
      return
      end
