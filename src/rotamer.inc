
c  This file holds the variable needed for rotamer library moves

c     the maximum number of rotamer dihedrals per residue
      integer MAXROTAMERDIHS
      parameter(MAXROTAMERDIHS=5)

c     the maximum number of rotamer states per residue
      integer MAXROTAMERSTATES
      parameter(MAXROTAMERSTATES=60)

c     the number of rotamer states for each type of residue
      integer NRotamerStates(MAXTEMPLATES)

c     the number of rotamer dihedrals for each type of residue
      integer NRotamerDihs(MAXTEMPLATES)

c     the names of the dihedral atoms
      character*4 RotamerDihName(MAXTEMPLATES,MAXROTAMERDIHS,4)

c     the actual values for each rotamer for each type of residue
      double precision RotamerVal(MAXTEMPLATES,MAXROTAMERSTATES,MAXROTAMERDIHS)

c     the number which identifies which rotamer library each residue in a 
c     protein uses. A value of 0 indicates that the residue does not have
c     rotamer library entry, or that it has an empty rotamer library entry
      integer ResRotamerID(MAXPROTEINS,MAXRESIDUES)

c     the identifying number that locates the rotamer dihedral for each
c     side chain atom in the residue
      integer RotamerID(MAXPROTEINS,MAXRESIDUES,MAXSCATOMS)

c     the probability of doing a big dihedral move instead of a rotamer move
c     (a 10% chance of a big dihedral move would be a value of 0.1)
      double precision BigMoveProb

c     the delta for big moves, typically equal to 180 degrees (for +/- 180 degree moves)
      double precision BigMoveDelta

c     the precision for rotamer moves, typically equal to 10 degrees (e.g. we move
c     the dihedral to be within +/- 10 degrees of the rotamer value)
      double precision RotamerDelta

      common /ROTAMER/ BigMoveProb, BigMoveDelta, RotamerDelta,
     .                 RotamerVal, NRotamerDihs, NRotamerStates,
     .                 ResRotamerID, RotamerID, RotamerDihName

      save /ROTAMER/
