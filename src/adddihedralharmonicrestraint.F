      subroutine addDihedralHarmonicRestraint(itype1,imol1,ires1,iatm1,
     .  itype2,imol2,ires2,iatm2,itype3,imol3,ires3,iatm3,itype4,
     .  imol4,ires4,iatm4,eq,k)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'printer.inc'
      include 'extraenergies.inc'
c#############################################
c     
c     Add a dihedral harmonic restraint ! 
c
c     Julien Michel February 2005
c
c#############################################
      integer itype1,imol1,iatm1,itype2,imol2,iatm2
      integer itype3,imol3,iatm3,itype4,imol4,iatm4
      integer ires1,ires2,ires3,ires4
      double precision eq,k
      
      nDihedralHarmonicRestraints = nDihedralHarmonicRestraints + 1
      DihedralHarmEq(nDihedralHarmonicRestraints) = eq
      DihedralHarmConst(nDihedralHarmonicRestraints) = k
      DihedralHarmAtom1(nDihedralHarmonicRestraints,1) = itype1
      DihedralHarmAtom1(nDihedralHarmonicRestraints,2) = imol1
      DihedralHarmAtom1(nDihedralHarmonicRestraints,3) = ires1
      DihedralHarmAtom1(nDihedralHarmonicRestraints,4) = iatm1

      if (itype2 .eq. 0) then
        ZmatDihedralHarm(nDihedralHarmonicRestraints) = .true.
c       Set other harm atom to 0 for 'safety'
        DihedralHarmAtom2(nDihedralHarmonicRestraints,1) = 0
        DihedralHarmAtom2(nDihedralHarmonicRestraints,2) = 0
        DihedralHarmAtom2(nDihedralHarmonicRestraints,3) = 0
        DihedralHarmAtom2(nDihedralHarmonicRestraints,4) = 0
        DihedralHarmAtom3(nDihedralHarmonicRestraints,1) = 0
        DihedralHarmAtom3(nDihedralHarmonicRestraints,2) = 0
        DihedralHarmAtom3(nDihedralHarmonicRestraints,3) = 0
        DihedralHarmAtom3(nDihedralHarmonicRestraints,4) = 0
        DihedralHarmAtom4(nDihedralHarmonicRestraints,1) = 0
        DihedralHarmAtom4(nDihedralHarmonicRestraints,2) = 0
        DihedralHarmAtom4(nDihedralHarmonicRestraints,3) = 0
        DihedralHarmAtom4(nDihedralHarmonicRestraints,4) = 0
      else
c       The torsion is defined by a zmat
        DihedralHarmAtom2(nDihedralHarmonicRestraints,1) = itype2
        DihedralHarmAtom2(nDihedralHarmonicRestraints,2) = imol2
        DihedralHarmAtom2(nDihedralHarmonicRestraints,3) = ires2
        DihedralHarmAtom2(nDihedralHarmonicRestraints,4) = iatm2
        DihedralHarmAtom3(nDihedralHarmonicRestraints,1) = itype3
        DihedralHarmAtom3(nDihedralHarmonicRestraints,2) = imol3
        DihedralHarmAtom3(nDihedralHarmonicRestraints,3) = ires3
        DihedralHarmAtom3(nDihedralHarmonicRestraints,4) = iatm3
        DihedralHarmAtom4(nDihedralHarmonicRestraints,1) = itype4
        DihedralHarmAtom4(nDihedralHarmonicRestraints,2) = imol4
        DihedralHarmAtom4(nDihedralHarmonicRestraints,3) = ires4
        DihedralHarmAtom4(nDihedralHarmonicRestraints,4) = iatm4
      endif
      

      return
      end
