c################################################
c
c     This file contains the dimensions 
c     for the program (e.g. maximum number
c     of proteins, MAXPROTEINS etc.)
c
c     This file *must* be present as the first
c     include file of any routine/function that
c     includes any other include files.
c
c     (C) Christopher Woods, October 2004
c
c###############################################

      implicit none

c parameters
c     maximum number of proteins to load
      integer MAXPROTEINS
      parameter(MAXPROTEINS=3)
c     the null string
      character*4 NULL
      parameter(NULL="NULL")
      character*4 NULLCLJ
      parameter(NULLCLJ="-?")
c     maximum number of residues per protein
      integer MAXRESIDUES
      parameter(MAXRESIDUES=1000)
c     maximum number of bonded residues per protein
      integer MAXRESRESBONDS
      parameter(MAXRESRESBONDS=MAXRESIDUES*2)
c     maximum number of sidechain atoms per residue
      integer MAXSCATOMS
      parameter(MAXSCATOMS=30)
c     the maximum possible atoms for each protein 
c     - should equal MAXRESIDUES*(MAXSCATOMS+3)
      integer MAXATOMS
      parameter(MAXATOMS=MAXRESIDUES*(MAXSCATOMS+3))
c     the maximum numbers of bonds, angles, dihedrals and urey-bradleys per atom
      integer MAXBONDSPERRESIDUE
      parameter(MAXBONDSPERRESIDUE=30)
      integer MAXANGLESPERRESIDUE
      parameter(MAXANGLESPERRESIDUE=100)
      integer MAXDIHEDRALSPERRESIDUE
      parameter(MAXDIHEDRALSPERRESIDUE=100)
      integer MAXUREYBRADLEYSPERRESIDUE
      parameter(MAXUREYBRADLEYSPERRESIDUE=100)
c     the total maximum number of bonds/angles/dihedrals in the protein
      integer MAXBONDS
      parameter(MAXBONDS=(MAXBONDSPERRESIDUE+5)*MAXRESIDUES)
      integer MAXANGLES
      parameter(MAXANGLES=(MAXANGLESPERRESIDUE+5)*MAXRESIDUES)
      integer MAXDIHEDRALS
      parameter(MAXDIHEDRALS=(MAXDIHEDRALSPERRESIDUE+5)*MAXRESIDUES)
      integer MAXUREYBRADLEYS
      parameter(MAXUREYBRADLEYS=(MAXUREYBRADLEYSPERRESIDUE+5)*MAXRESIDUES)
c     the maximum number of inter- and intra residue nb and 14 pairs
      integer MAXRESNB
      parameter(MAXRESNB=MAXSCATOMS*7)
      integer MAXRESNB14
      parameter(MAXRESNB14=MAXSCATOMS*2)
      integer MAXRESRESNB
      parameter(MAXRESRESNB=MAXSCATOMS)
      integer MAXRESRESNB14
      parameter(MAXRESRESNB14=MAXSCATOMS)
c     the maximum number of acidic residues per protein
      integer MAXACIDICRES
      parameter(MAXACIDICRES=20)
c     the maximum number of basic residues per protein
      integer MAXBASICRES
      parameter(MAXBASICRES=20)
c     the total maximum number of charge/LJ parameters
      integer MAXCLJ
      parameter(MAXCLJ=10000)
c     the maximum number of different CLJ parameters within a single simulation
      integer MAXCLJTYPES
      parameter(MAXCLJTYPES=500)
c     the maximum number of variable parameters per simulation...
      integer MAXVARIABLEPARAM
      parameter(MAXVARIABLEPARAM=100)
c     the total maximum number of pairs of different charge/LJ parameters 
c     in a single simulation. (should equal MAXCLJTYPES**2 / 2)
      integer MAXCLJPAIRS
      parameter(MAXCLJPAIRS=(MAXCLJTYPES**2)/2)
c     the maximum number of bond types
      integer MAXBNDPARAM
      parameter(MAXBNDPARAM=5000)
c     the maximum number of angle types
      integer MAXANGPARAM
      parameter(MAXANGPARAM=5100)
c     the maximum number of dihedral types
      integer MAXDIHPARAM
      parameter(MAXDIHPARAM=5000)
c     the maximum number of dihedral terms
      integer MAXDIHTERMS
      parameter(MAXDIHTERMS=6000)
c     the maximum number of terms in each dihedral
      integer MAXDIHTERMSPERDIHEDRAL
      parameter(MAXDIHTERMSPERDIHEDRAL=6)
c     the maximum number of urey bradley types
      integer MAXUBYPARAM
      parameter(MAXUBYPARAM=500)

c     now the parameters corresponding to the solvents...
c     first the maximum number of solvents
      integer MAXSOLVENTS
      parameter(MAXSOLVENTS=75000)
c     now the maximum number of atoms per solvent molecule
      integer MAXSOLVENTATOMS
      parameter(MAXSOLVENTATOMS=10)

c     the maximum number of templates, with maximum numbers of 
c     bonds/angle/dihedrals in each template...
      integer MAXTEMPLATES
      parameter(MAXTEMPLATES=500)
      integer MAXTEMPLATEBONDS
      parameter(MAXTEMPLATEBONDS=100)
      integer MAXTEMPLATEANGLES
      parameter(MAXTEMPLATEANGLES=100)
      integer MAXTEMPLATEDIHEDRALS
      parameter(MAXTEMPLATEDIHEDRALS=100)
      integer MAXTEMPLATEUREYBRADLEYS
      parameter(MAXTEMPLATEUREYBRADLEYS=100)
      integer MAXTEMPLATEATMS
      parameter(MAXTEMPLATEATMS=100)
      integer MAXCHAINTEMPLATES
      parameter(MAXCHAINTEMPLATES=100)
      integer MAXCHAINATOMS
      parameter(MAXCHAINATOMS=10)
      integer MAXCHAINBONDS
      parameter(MAXCHAINBONDS=10)
      integer MAXCHAINANGLES
      parameter(MAXCHAINANGLES=10)
      integer MAXCHAINDIHEDRALS
      parameter(MAXCHAINDIHEDRALS=10)
      integer MAXCHAINUREYBRADLEYS
      parameter(MAXCHAINUREYBRADLEYS=10)
c     the maximum number of override templates
      integer MAXOVERRIDES
      parameter(MAXOVERRIDES=MAXTEMPLATES*4)
      integer MAXOVERRIDEPARAMS
      parameter(MAXOVERRIDEPARAMS=MAXTEMPLATEATMS)
c     the same for the molecule templates
      integer MAXMOLTEMPLATES
      parameter(MAXMOLTEMPLATES=150)
      integer MAXMOLTEMPLATEBONDS
      parameter(MAXMOLTEMPLATEBONDS=200)
      integer MAXMOLTEMPLATEANGLES
      parameter(MAXMOLTEMPLATEANGLES=200)
      integer MAXMOLTEMPLATEDIHEDRALS
      parameter(MAXMOLTEMPLATEDIHEDRALS=200)
      integer MAXMOLTEMPLATEUREYBRADLEYS
      parameter(MAXMOLTEMPLATEUREYBRADLEYS=200)
      integer MAXMOLTEMPLATEATMS
      parameter(MAXMOLTEMPLATEATMS=200)
      integer MAXMOLTEMPLVAR
      parameter(MAXMOLTEMPLVAR=MAXTEMPLATEATMS)
c     now the same for the solvent templates
      integer MAXSVNTEMPLATES
      parameter(MAXSVNTEMPLATES=150)
      integer MAXSVNTEMPLATEATMS
      parameter(MAXSVNTEMPLATEATMS=MAXSOLVENTATOMS)
c     the value of a NULL name
      character*4 NULLNAME
      parameter(NULLNAME = "    ")
c     now the parameters corresponding to the solutes...
      integer MAXSOLUTES
      parameter(MAXSOLUTES=60)
      integer MAXSOLUTERESIDUES
      parameter(MAXSOLUTERESIDUES=10)
      integer MAXSOLUTEATOMSPERRESIDUE
      parameter(MAXSOLUTEATOMSPERRESIDUE=100)
      integer MAXSOLUTEATOMS
      parameter(MAXSOLUTEATOMS=MAXSOLUTERESIDUES*MAXSOLUTEATOMSPERRESIDUE)
      integer MAXSOLVARIABLES
      parameter(MAXSOLVARIABLES=MAXSOLUTEATOMS)
      integer MAXSOLUTEBONDS
      parameter(MAXSOLUTEBONDS=500)
      integer MAXSOLUTEANGLES
      parameter(MAXSOLUTEANGLES=1000)
      integer MAXSOLUTEDIHEDRALS
      parameter(MAXSOLUTEDIHEDRALS=2000)
      integer MAXSOLUTEUREYBRADLEYS
      parameter(MAXSOLUTEUREYBRADLEYS=1000)
      integer MAXSOLUTENBPAIRS
      parameter(MAXSOLUTENBPAIRS=10000)
      
c     now the parameters relating to any stacks used in the code
      integer MAXRESSTACKDEPTH
      parameter(MAXRESSTACKDEPTH=30)
      integer MAXRESCOORDDEPTH
      parameter(MAXRESCOORDDEPTH=4)
      integer MAXSOLCOORDDEPTH
      parameter(MAXSOLCOORDDEPTH=4)
      integer MAXSVNCOORDDEPTH
      parameter(MAXSVNCOORDDEPTH=4)
      integer MAXSVNSTACKDEPTH
      parameter(MAXSVNSTACKDEPTH=10000000)
      integer MAXCAVCOORDEPTH
      parameter(MAXCAVCOORDEPTH=1000)

c     grid bits
      integer MAXXGRID
      parameter(MAXXGRID=10000000)
      integer MAXYGRID
      parameter(MAXYGRID=10000000)
      integer MAXZGRID
      parameter(MAXZGRID=10000000)
      integer MAXRESERVOIR
      parameter(MAXRESERVOIR=50000)
      
c     Parallel tempering bits
c     Maximum allowed number of lambdas
      integer MAXLAMBDAS
      parameter (MAXLAMBDAS=50)
c     Maximum allowed number of temperatures
      integer MAXTEMPS
      parameter (MAXTEMPS=50) 
c     Maximum allowed number of adams values
      integer MAXADAMS
      parameter (MAXADAMS=50)
c     Maximum allowed number of jaws2 bias
      integer MAXJAWS2
      parameter (MAXJAWS2=50)
c     Maximum allowed number of ladders lambdas/temperatures/adams/jaws2
      integer MAXLADDERS
      parameter (MAXLADDERS=50)

c     Maximum number of possible dump statements
      integer MAXDUMP
      parameter (MAXDUMP=25)
