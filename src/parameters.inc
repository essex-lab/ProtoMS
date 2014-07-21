
c  This file holds the variables that make up the parameter database

c characters
c     the amber names of the atoms in the parameter database...
      character*2 AtmDB(MAXCLJ)
c     the amber names of the variable parameters
      character*2 VariableAtm(MAXVARIABLEPARAM,2)
c     the names of the bond parameters (AMBER name)
      character*2 BndName(MAXBNDPARAM,2)
c     the names of the angle parameters (AMBER name)
      character*2 AngName(MAXANGPARAM,3)
c     the names of the dihedral parameters (AMBER name)
      character*2 DihName(MAXDIHPARAM,4)
c     the names of the urey-bradley parameters (AMBER name)
      character*2 UByName(MAXUBYPARAM,3)
            
      common /CPARAMDB/ VariableAtm,AtmDB,BndName,AngName,DihName,UByName
      save /CPARAMDB/

c     whether or not this clj type is loaded in the parameter database, and 
c     whether it is being used in the simulations
      logical GotCLJ(MAXCLJ,2)
c     whether or not the bond/angle/dihedral/uby parameters have been loaded
      logical GotBndParam(MAXBNDPARAM)
      logical GotAngParam(MAXANGPARAM)
      logical GotDihParam(MAXDIHPARAM)
      logical GotUByParam(MAXUBYPARAM)
c     whether or not this dihedral term is loaded in the parameter database
      logical GotDihTerm(MAXDIHTERMS)

      common /LPARAMDB/ GotCLJ,GotBndParam,GotAngParam,GotDihParam,
     .          GotUByParam,GotDihTerm
      
      save /LPARAMDB/
      
c     the list of clj types used in the simulation...
      integer CLJType(MAXCLJTYPES)
c     the index in the clj pair array, for pair of parameters i,j
      integer CLJPair(MAXCLJTYPES,MAXCLJTYPES)
c     the random number seed used to generate the trajectory
      integer ranseed
c     the number of variable atoms in the simulation
      integer NVariableCLJ
      integer NVarBondParam
      integer NVarAngleParam
      integer NVarDihedralParam
      integer NVarUreyBradleyParam
c     the proton number of each clj parameter
      integer ProtonNumDB(MAXCLJ)
c     the proton number of variable parameters
      integer VariableProtNum(MAXVARIABLEPARAM,2)
c     the list of variable clj parameters used in the simulation
c     (the parameters that vary with lambda)
      integer VariableCLJ(MAXVARIABLEPARAM,2)
c     the index in the bond parameter array for bond name i
      integer BndParamIndex(MAXBNDPARAM)
c     the index in the angle parameter array for angle name i
      integer AngParamIndex(MAXANGPARAM)
c     the index in the dihedral parameter array for dihedral name i
      integer DihParamIndex(MAXDIHPARAM)
c     the number of dihedral terms for this indexed dihedral
      integer NDihTerms(MAXDIHPARAM)
c     the list of indexes into the dihedral parameter array for dihedral i
      integer DihTermList(MAXDIHPARAM,MAXDIHTERMSPERDIHEDRAL)
c     parameters for variable dihedrals - the indicies for lambda=0 and lambda=1
      integer VariableDihParam(MAXVARIABLEPARAM,2)
c     the index in the urey-bradley parameter array for urey-bradley name i
      integer UByParamIndex(MAXUBYPARAM)
c     the number of cljtypes and pairs in the simulation...
      integer NCLJTypes
c     the number of bond and angle parameters
      integer NBndParam,NAngParam,NDihParam,NUByParam

      common /IPARAMDB/ CLJType,NCLJTypes,NBndParam,NAngParam,CLJPair,
     .          NDihParam,NUByParam,VariableCLJ,VariableProtNum,NVariableCLJ,
     .          VariableDihParam,ProtonNumDB,DihParamIndex,ranseed,
     .          NVarBondParam,NVarAngleParam,NVarDihedralParam,NVarUreyBradleyParam,
     .          BndParamIndex,AngParamIndex,NDihTerms,DihTermList,UByParamIndex
      
      save /IPARAMDB/
      
c     the charge and LJ parameters loaded into the parameter database
      double precision ChgDB(MAXCLJ),SigDB(MAXCLJ),EpsDB(MAXCLJ)
c     multiples of pairs of charge and LJ parameters
      double precision ChgIJ(MAXCLJPAIRS),AIJ(MAXCLJPAIRS)
      double precision BIJ(MAXCLJPAIRS)

c     the variable charge and LJ parameters, at lam=0, lam=lambda
c     and lam=1.0
      double precision VariableChg(MAXVARIABLEPARAM,5)
      double precision VariableSig(MAXVARIABLEPARAM,5)
      double precision VariableEps(MAXVARIABLEPARAM,5)
c     bond parameters
      double precision BndParam(MAXBNDPARAM,2)
c     angle parameters
      double precision AngParam(MAXANGPARAM,2)
c     urey-bradley parameters
      double precision UByParam(MAXUBYPARAM,2)
c     dihedral parameters
      double precision DihTerm(MAXDIHPARAM,4)
c     the variable bond/angle/dihedral parameters
      double precision VariableBndParam(MAXVARIABLEPARAM,5,2)
      double precision VariableAngParam(MAXVARIABLEPARAM,5,2)
      double precision VariableUByParam(MAXVARIABLEPARAM,5,2)

      common /DPARAMDB/ ChgDB,SigDB,EpsDB,ChgIJ,AIJ,BIJ,DihTerm,
     .           BndParam,AngParam,UByParam,
     .           VariableBndParam,VariableAngParam,VariableUByParam,
     .           VariableChg,VariableSig,VariableEps
     
      save /DPARAMDB/
