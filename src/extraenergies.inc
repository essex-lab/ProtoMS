
c   This file holds the variables that are associated with the extra energy terms
c   Place variables that you need for the extra energies here (or you can use your
c   own include file)

c     the total extra energy variables
      double precision oldExtraEnergy,avgExtraEnergy,avgExtraEnergy2,
     .                 oldExtraEnergyB,avgExtraEnergyB,avgExtraEnergyB2,
     .                 oldExtraEnergyF,avgExtraEnergyF,avgExtraEnergyF2,
     .                 oldExtraEnergyPart,newExtraEnergyPart,
     .                 oldExtraEnergyPartB,newExtraEnergyPartB,
     .                 oldExtraEnergyPartF,newExtraEnergyPartF
     
      common /EXTRAENERGIES/ oldExtraEnergy,avgExtraEnergy,avgExtraEnergy2,
     .                 oldExtraEnergyB,avgExtraEnergyB,avgExtraEnergyB2,
     .                 oldExtraEnergyF,avgExtraEnergyF,avgExtraEnergyF2,
     .                 oldExtraEnergyPart,newExtraEnergyPart,
     .                 oldExtraEnergyPartB,newExtraEnergyPartB,
     .                 oldExtraEnergyPartF,newExtraEnergyPartF
     
c     variables associated with the solvent cap
      double precision oldSvnCapEnergy,avgSvnCapEnergy,avgSvnCapEnergy2,
     .                 oldSvnCapEnergyPart,newSvnCapEnergyPart
    
      common /EXTRACAP/ oldSvnCapEnergy,avgSvnCapEnergy,avgSvnCapEnergy2,
     .                  oldSvnCapEnergyPart,newSvnCapEnergyPart
    
      save /EXTRACAP/

c     variables associated with harmonic restraints
      integer NHarmonicRestraints
      integer MAXHARMONICRESTRAINTS
      parameter(MAXHARMONICRESTRAINTS=100)

c     the atom codes for each of the harmonic restraints - contains the
c     molecule type, molecule number and atom number for each restraints
c     (see getCoordinates() for more info)
      integer HarmAtom(MAXHARMONICRESTRAINTS,3)      

c     the point in space that the atom is restrained to
      double precision HarmPoint(MAXHARMONICRESTRAINTS,3)      
                  
c     the harmonic force constant
      double precision HarmConst(MAXHARMONICRESTRAINTS)      
      
c     whether to scale a restraint by lambda
      logical HarmScale(MAXHARMONICRESTRAINTS)

      double precision oldHarmEnergy,avgHarmEnergy,avgHarmEnergy2,
     .                 oldHarmEnergyB,avgHarmEnergyB,avgHarmEnergyB2,
     .                 oldHarmEnergyF,avgHarmEnergyF,avgHarmEnergyF2,
     .                 oldHarmEnergyPart,newHarmEnergyPart,
     .                 oldHarmEnergyPartB,newHarmEnergyPartB,
     .                 oldHarmEnergyPartF,newHarmEnergyPartF
     
      common /EXTRAHARM/ oldHarmEnergy,avgHarmEnergy,avgHarmEnergy2,
     .                 oldHarmEnergyB,avgHarmEnergyB,avgHarmEnergyB2,
     .                 oldHarmEnergyF,avgHarmEnergyF,avgHarmEnergyF2,
     .                 oldHarmEnergyPart,newHarmEnergyPart,
     .                 oldHarmEnergyPartB,newHarmEnergyPartB,
     .                 oldHarmEnergyPartF,newHarmEnergyPartF,
     .                 HarmPoint,HarmConst,HarmAtom,NHarmonicRestraints,
     .                 HarmScale


c    variables associated with bond harmonic restraints
      integer NBndHarmonicRestraints
      integer MAXBNDHARMONICRESTRAINTS
      parameter(MAXBNDHARMONICRESTRAINTS=300)
      integer BndHarmAtom1(MAXBNDHARMONICRESTRAINTS,3)         
      integer BndHarmAtom2(MAXBNDHARMONICRESTRAINTS,3)         
      double precision BndHarmConst(MAXBNDHARMONICRESTRAINTS)      
      double precision BndHarmEq(MAXBNDHARMONICRESTRAINTS)      

      double precision oldBndHarmEnergy,avgBndHarmEnergy,avgBndHarmEnergy2,
     .                 oldBndHarmEnergyB,avgBndHarmEnergyB,avgBndHarmEnergyB2,
     .                 oldBndHarmEnergyF,avgBndHarmEnergyF,avgBndHarmEnergyF2,
     .                 oldBndHarmEnergyPart,newBndHarmEnergyPart,
     .                 oldBndHarmEnergyPartB,newBndHarmEnergyPartB,
     .                 oldBndHarmEnergyPartF,newBndHarmEnergyPartF

      common /EXTRABNDHARM/ oldBndHarmEnergy,avgBndHarmEnergy,avgBndHarmEnergy2,
     .                 oldBndHarmEnergyB,avgBndHarmEnergyB,avgBndHarmEnergyB2,
     .                 oldBndHarmEnergyF,avgBndHarmEnergyF,avgBndHarmEnergyF2,
     .                 oldBndHarmEnergyPart,newBndHarmEnergyPart,
     .                 oldBndHarmEnergyPartB,newBndHarmEnergyPartB,
     .                 oldBndHarmEnergyPartF,newBndHarmEnergyPartF,
     .                 BndHarmConst,BndHarmEq,BndHarmAtom1,BndHarmAtom2,NBndHarmonicRestraints

c    variables associated with torsional harmonic restraints
      integer NDihedralHarmonicRestraints
      integer MAXDIHEDRALHARMONICRESTRAINTS
      parameter(MAXDIHEDRALHARMONICRESTRAINTS=10)
c    The four atoms involved in the dihedral
c     what goes in is ityp,imol,ires,iatm
      integer DihedralHarmAtom1(MAXDIHEDRALHARMONICRESTRAINTS,4)
      integer DihedralHarmAtom2(MAXDIHEDRALHARMONICRESTRAINTS,4)
      integer DihedralHarmAtom3(MAXDIHEDRALHARMONICRESTRAINTS,4)
      integer DihedralHarmAtom4(MAXDIHEDRALHARMONICRESTRAINTS,4)
c    Whether or not the dihedral is actually the one defined in a zmatrix
      logical ZmatDihedralHarm(MAXDIHEDRALHARMONICRESTRAINTS)
c    The value the dihedral is restrained to
      double precision DihedralHarmEq(MAXDIHEDRALHARMONICRESTRAINTS)
c    The force constant used to hold it in place
      double precision DihedralHarmConst(MAXDIHEDRALHARMONICRESTRAINTS)  
   
      double precision oldDihHarmEnergy,avgDihHarmEnergy,avgDihHarmEnergy2,
     .                 oldDihHarmEnergyB,avgDihHarmEnergyB,avgDihHarmEnergyB2,
     .                 oldDihHarmEnergyF,avgDihHarmEnergyF,avgDihHarmEnergyF2,
     .                 newDihHarmEnergyPart,newDihHarmEnergyPartF,
     .                   newDihHarmEnergyPartB,oldDihHarmEnergyPart,
     .                   oldDihHarmEnergyPartB,oldDihHarmEnergyPartF 

      common /EXTRADIHHARM/oldDihHarmEnergy,avgDihHarmEnergy,avgDihHarmEnergy2,
     .                 oldDihHarmEnergyB,avgDihHarmEnergyB,avgDihHarmEnergyB2,
     .                 oldDihHarmEnergyF,avgDihHarmEnergyF,avgDihHarmEnergyF2,
     .                 newDihHarmEnergyPart,newDihHarmEnergyPartF,
     .                   newDihHarmEnergyPartB,oldDihHarmEnergyPart,
     .                   oldDihHarmEnergyPartB,oldDihHarmEnergyPartF,
     .                   DihedralHarmEq,DihedralHarmConst,
     .                   NDihedralHarmonicRestraints,DihedralHarmAtom1,
     .                   DihedralHarmAtom2,DihedralHarmAtom3,DihedralHarmAtom4,
     .                   ZmatDihedralHarm
c      Cartesian Harmonic Flat Bottom Restraints
c     variables associated with harmonic restraints
      integer NFlatBottomRestraints
      integer MAXFLATBOTTOMRESTRAINTS
      parameter(MAXFLATBOTTOMRESTRAINTS=30)

c     the atom codes for each of the harmonic restraints - contains the
c     molecule type, molecule number and atom number for each restraints
c     (see getCoordinates() for more info)
      integer FBotAtom(MAXFLATBOTTOMRESTRAINTS,3)      

c     the point in space that the atom is restrained to
      double precision FBotPoint(MAXFLATBOTTOMRESTRAINTS,3)      
                  
c     the harmonic force constant
      double precision FBotConst(MAXFLATBOTTOMRESTRAINTS)      
      double precision FBotWidth(MAXFLATBOTTOMRESTRAINTS)

      double precision oldFBotEnergy,avgFBotEnergy,avgFBotEnergy2,
     .                 oldFBotEnergyB,avgFBotEnergyB,avgFBotEnergyB2,
     .                 oldFBotEnergyF,avgFBotEnergyF,avgFBotEnergyF2,
     .                 oldFBotEnergyPart,newFBotEnergyPart,
     .                 oldFBotEnergyPartB,newFBotEnergyPartB,
     .                 oldFBotEnergyPartF,newFBotEnergyPartF
     
      common /EXTRAFBOT/ oldFBotEnergy,avgFBotEnergy,avgFBotEnergy2,
     .                 oldFBotEnergyB,avgFBotEnergyB,avgFBotEnergyB2,
     .                 oldFBotEnergyF,avgFBotEnergyF,avgFBotEnergyF2,
     .                 oldFBotEnergyPart,newFBotEnergyPart,
     .                 oldFBotEnergyPartB,newFBotEnergyPartB,
     .                 oldFBotEnergyPartF,newFBotEnergyPartF,
     .                   FBotWidth,FBotPoint,FBotConst,FBotAtom,
     .                   NFlatBottomRestraints

c      dihedral LamTor restraints
      integer nDihedralLamTorRestraints
      integer MAXDIHEDRALLAMTORRESTRAINTS
      parameter(MAXDIHEDRALLAMTORRESTRAINTS=10)
c    The four atoms involved in the dihedral
c     what goes in is ityp,imol,ires,iatm
      integer DihedralLamTorAtom1(MAXDIHEDRALLAMTORRESTRAINTS,4)
      integer DihedralLamTorAtom2(MAXDIHEDRALLAMTORRESTRAINTS,4)
      integer DihedralLamTorAtom3(MAXDIHEDRALLAMTORRESTRAINTS,4)
      integer DihedralLamTorAtom4(MAXDIHEDRALLAMTORRESTRAINTS,4)
c    Whether or not the dihedral is actually the one defined in a zmatrix
      logical ZmatDihedralLamTor(MAXDIHEDRALLAMTORRESTRAINTS)
c    The value the dihedral is restrained to
      double precision DihedralLamTorEq(MAXDIHEDRALLAMTORRESTRAINTS)
      double precision DihedralLamTorConst(MAXDIHEDRALLAMTORRESTRAINTS)  
      double precision DihedralLamTorNforw(MAXDIHEDRALLAMTORRESTRAINTS)  
      double precision DihedralLamTorNback(MAXDIHEDRALLAMTORRESTRAINTS)          

      double precision oldDihLamTorEnergy,avgDihLamTorEnergy,avgDihLamTorEnergy2,
     .                 oldDihLamTorEnergyB,avgDihLamTorEnergyB,avgDihLamTorEnergyB2,
     .                 oldDihLamTorEnergyF,avgDihLamTorEnergyF,avgDihLamTorEnergyF2,
     .                 newDihLamTorEnergyPart,newDihLamTorEnergyPartF,
     .                   newDihLamTorEnergyPartB,oldDihLamTorEnergyPart,
     .                   oldDihLamTorEnergyPartB,oldDihLamTorEnergyPartF 

      common /EXTRADIHLAMTOR/oldDihLamTorEnergy,avgDihLamTorEnergy,avgDihLamTorEnergy2,
     .                 oldDihLamTorEnergyB,avgDihLamTorEnergyB,avgDihLamTorEnergyB2,
     .                 oldDihLamTorEnergyF,avgDihLamTorEnergyF,avgDihLamTorEnergyF2,
     .                 newDihLamTorEnergyPart,newDihLamTorEnergyPartF,
     .                   newDihLamTorEnergyPartB,oldDihLamTorEnergyPart,
     .                   oldDihLamTorEnergyPartB,oldDihLamTorEnergyPartF,
     .                   DihedralLamTorEq,DihedralLamTorConst,
     .                 DihedralLamTorNforw,DihedralLamTorNback,
     .                   NDihedralLamTorRestraints,DihedralLamTorAtom1,
     .                   DihedralLamTorAtom2,DihedralLamTorAtom3,DihedralLamTorAtom4,
     .                   ZmatDihedralLamTor
c      GBSA extension
      include 'surface.inc'
      include 'gb.inc'

      integer NHardwallPotential
      integer MAXHARDWALLATOMS
      parameter(MAXHARDWALLATOMS=100)

c     the id for the molecule of interest
      integer HardwallMolecule(MAXHARDWALLATOMS)

c     the point in space that the atom is restrained to
      double precision HardwallCenter(MAXHARDWALLATOMS,3)

c     the cutoff radius
      double precision HardwallRadius(MAXHARDWALLATOMS)

c     the energy added if within the radius
      double precision HardwallextraEnergy(MAXHARDWALLATOMS)


      double precision oldHardwallEnergy,avgHardwallEnergy,avgHardwallEnergy2,
     .                 oldHardwallEnergyB,avgHardwallEnergyB,avgHardwallEnergyB2,
     .                 oldHardwallEnergyF,avgHardwallEnergyF,avgHardwallEnergyF2,
     .                 oldHardwallEnergyPart,newHardwallEnergyPart,
     .                 oldHardwallEnergyPartB,newHardwallEnergyPartB,
     .                 oldHardwallEnergyPartF,newHardwallEnergyPartF

      common /EXTRAHARDWALL/ oldHardwallEnergy,avgHardwallEnergy,avgHardwallEnergy2,
     .                 oldHardwallEnergyB,avgHardwallEnergyB,avgHardwallEnergyB2,
     .                 oldHardwallEnergyF,avgHardwallEnergyF,avgHardwallEnergyF2,
     .                 oldHardwallEnergyPart,newHardwallEnergyPart,
     .                 oldHardwallEnergyPartB,newHardwallEnergyPartB,
     .                 oldHardwallEnergyPartF,newHardwallEnergyPartF,
     .                 HardwallCenter,HardwallRadius,HardwallMolecule,HardwallextraEnergy,NHardwallPotential

      
