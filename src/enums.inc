
c##########################################
c
c     These are the 'enums' or '#define's
c     used in ProtoMC.
c
c     (C) Christopher Woods, October 2004
c
c##########################################

c     type of monte carlo test
      integer NPTMCTEST
      integer GCMCTESTA
      integer GCMCTESTD
      parameter(NPTMCTEST=1)
      parameter (GCMCTESTA=2)
      parameter (GCMCTESTD=3)
                                     
c     a bond/angle/dihedral
      integer BONDS,ANGLES,DIHEDRALS,UREYBRADLEYS
      parameter(BONDS=1)
      parameter(ANGLES=2)
      parameter(DIHEDRALS=3)
      parameter(UREYBRADLEYS=4)

c     position in a chain, first, middle, last, single or all positions
      integer NOPOS,FIRST,MIDDLE,LAST,SINGLE,ALLPOS
      parameter(NOPOS=0)
      parameter(FIRST=1)
      parameter(MIDDLE=2)
      parameter(LAST=3)
      parameter(SINGLE=4)
      parameter(ALLPOS=5)
            
c     now the different parts of the system
      integer PROTEIN,SOLUTE,SOLVENT,RESIDUE,ATOMS,GCSOLUTE
      parameter(PROTEIN = 1)
      parameter(SOLUTE = 2)
      parameter(SOLVENT = 3)
      parameter(RESIDUE = 4)
      parameter(ATOMS = 5)
      parameter(GCSOLUTE = 6)
c     old vs. new vs. don't care
      integer OLD,NEW,DONTCARE
      parameter(OLD = 1)
      parameter(NEW = 2)
      parameter(DONTCARE = 3)
c     integer key used to update the energies
      integer UPDATEENERGY
      parameter(UPDATEENERGY = 1)
      integer NOUPDATE
      parameter(NOUPDATE = 0)
c     now the different types of move
      integer TRYSOLUTE,TRYSOLVENT,TRYPROTEIN,TRYTITRATE,TRYTHETA,TRYVOLUME,TRYROTAMER
      integer TRYINSERTION,TRYDELETION,TRYGCSOLUTE,TRYSAMPLE
      parameter(TRYSOLUTE=1)
      parameter(TRYSOLVENT=2)
      parameter(TRYPROTEIN=3)
      parameter(TRYTITRATE=4)
      parameter(TRYTHETA=5)
      parameter(TRYVOLUME=6)
      parameter(TRYROTAMER=7)
      parameter(TRYINSERTION=8)
      parameter(TRYDELETION=9)
      parameter(TRYGCSOLUTE=10)
      parameter(TRYSAMPLE=11)
c     now the different types of extra energy
      integer EXTRA_UMBRELLA,EXTRA_POSITION_RESTRAINT,EXTRA_BOND_RESTRAINT
      integer EXTRA_ANGLE_RESTRAINT,EXTRA_DIHEDRAL_RESTRAINT
      parameter(EXTRA_UMBRELLA=1)
      parameter(EXTRA_POSITION_RESTRAINT=2)
      parameter(EXTRA_BOND_RESTRAINT=3)
      parameter(EXTRA_ANGLE_RESTRAINT=4)
      parameter(EXTRA_DIHEDRAL_RESTRAINT=5)
      
c switches... (*) = default
c   the nonbonded cut test  (1)Protein-residue/Solute-molecule/solvent-molecule (*)
c                           (2)Protein-molecule/solute-molecule/solvent-molecule
      integer RESIDUEMOLECULE
      integer MOLECULEMOLECULE
      parameter(RESIDUEMOLECULE=1)
      parameter(MOLECULEMOLECULE=2)
c   the boundary type       (0)No boundaries
c                           (1)periodic boundaries (*)
c                           (2)Cubic harmonic boundaries (keep molecules in a box using harmonic potential)
c                           (3)Spherical harmonic boundaries
c
c                           (-1) Agree with whatever the solvent wants
      integer NOBOUNDARIES
      integer PERIODICBOUNDARIES
      integer CUBICHARMONIC
      integer SPHERICALHARMONIC
      integer AGREESOLVENT
      parameter(NOBOUNDARIES=0)
      parameter(PERIODICBOUNDARIES=1)
      parameter(CUBICHARMONIC=2)
      parameter(SPHERICALHARMONIC=3)
      parameter(AGREESOLVENT=-1)
       
