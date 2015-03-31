      subroutine getOptions
      include 'dimensions.inc'
      include 'metropolis.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'splitter.inc'
      include 'simparams.inc'
      include 'solvents.inc'
      include 'solutes.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'parameters.inc'
      include 'pressure.inc'
      include 'surface.inc'
      include 'gb.inc'
      include 'keyset.inc'
      include 'parallel.inc'
c#######################################################
c
c     This subroutine gets the values of most of the
c     user control options. The options are sets of
c     key-value pairs.
c
c     The values are set either in a command file
c     (see parsefile.F) or via environmental 
c     variables (see subroutine getenv2).
c
c     Christopher Woods, October 2004
c
c#######################################################

      integer i,j,lgth,lgth0,lgth1,lineLength,nkind,nold
      integer lowgroup,highgroup,k
      integer imol
c     the current system time and date            
      integer year,month,day,hour,mins,secs
      
      logical getLogical,getOption
      double precision getDouble
      integer getInt
      character*300 tempstring
      character*300 filename

      double precision delb,delf

      write(6,*)"GetOptions_start ",temperatureC


c     set initial, default values for options
      call setDefaultOptions

      write(6,*)"GetOptions_aftlambda ",temperatureC
            
c     Check to see if we should are debugging - this is the 'debug' option
      if (getOption('debug')) then
        protodebug = getLogical(words(1),protodebug)
        if (protodebug) then
          call printLine(INFO,"Printing debugging output")
C         turn on the debugging stream - initially to STDOUT
          streamhandle(DEBUG) = STDOUT
        endif
      endif

c     see if we should be pretty printing
      if (getOption('prettyprint')) prettyprint = getLogical(words(1),prettyprint)
            
c     check to see whether we will test the energy routines
      if (getOption('testenergy')) then
        testenergy = getLogical(words(1),testenergy)
        if (testenergy) call printLine(INFO,"Testing energy routines!")
      endif

c     get the current value of lambda
      if (getOption('lambda')) then
c       the format is 'lambda lambdaf lambdab'
        lambda = getDouble(words(1),lambda)
        if (nwords.ge.3) then
          lambdaf = getDouble(words(2),lambdaf)
          lambdab = getDouble(words(3),lambdab)
        endif      
      endif
      if (doLambdaRE) then
        lambda = lambdaval(mylambda)
        if (getOption('dlambda')) then
          deltalam = getDouble(words(1),0.0)
          lambdaf = lambda+deltalam
          lambdab = lambda-deltalam
          if (lambdaf.gt.1.0) lambdaf = 1.0
          if (lambdaf.lt.0.0) lambdaf = 0.0
          if (lambdab.gt.1.0) lambdab = 1.0
          if (lambdab.lt.0.0) lambdab = 0.0          
        endif
      endif

c     are we trying to match Sire?
      if (getOption('sire_compatibility_mode')) then
        SIRE_COMPATIBILITY_MODE = getLogical(words(1), .true.)
      
        if (SIRE_COMPATIBILITY_MODE) then
          call prettyLine(INFO,"Entering Sire compatibility mode")
        endif

      else
        SIRE_COMPATIBILITY_MODE = .false.
      endif

c     whether we are going to use chunks from protein pdb
      PdbParams = .false.
      if (getOption('pdbparams')) then
        PdbParams = getLogical(words(1), PdbParams)
        write(6,*) 'PdbParams'
        write(6,*) PdbParams
      endif  

c     whether jaws2 information on each move will be printed to results
      everyjaws2 = .false.
      if (getOption('everyjaws2')) then
        everyjaws2 = getLogical(words(1), everyjaws2)
        write(6,*) 'everyjaws2'
        write(6,*) everyjaws2
      endif    


c     whether free energy data will be printed in results
      if (getOption('printfe')) then
        if (words(1).eq.'bar') then
          printFE = .true.
          mbar = .false.
        else if (words(1).eq.'mbar') then
          printFE = .true.
          mbar = .true.
        endif
      else 
        printFE = .false.
      endif  

c     now get the type of boundaries to use...
      if (getOption('boundary')) then
c       the options are 'periodic' 'cubic' 'cap', 'none' or 'solvent'
        call lowercase(words(1))
        if (words(1).eq.'periodic') then
          boundary = PERIODICBOUNDARIES
          call printLine(INFO,"Requested periodic boundaries")
          if (nwords.ge.7) then
c           we can read the box dimensions from the option
            boxOrigin(1) = getDouble(words(2),ZERO)
            boxOrigin(2) = getDouble(words(3),ZERO)
            boxOrigin(3) = getDouble(words(4),ZERO)
            boxTop(1) = getDouble(words(5),ZERO)
            boxTop(2) = getDouble(words(6),ZERO)
            boxTop(3) = getDouble(words(7),ZERO)
            SolventBox = .true.
            SolventCap = .false.
          else if (nwords.ge.4) then
c           assume that only the box dimensions have been given, and that the box
c           is centered on the origin
            boxOrigin(1) = getDouble(words(2),ZERO)
            boxOrigin(2) = getDouble(words(3),ZERO)
            boxOrigin(3) = getDouble(words(4),ZERO)
            do i=1,3
              boxTop(i) = HALF*boxOrigin(i)
              boxOrigin(i) = -boxTop(i)
            enddo
          else
            call printLine(WARNING,"Estimating solvent box size. This will need to be equilibrated")
          endif                                  
        else if (words(1).eq.'cubic') then
          boundary = CUBICHARMONIC
          call printLine(INFO,
     .     "Requested cubic harmonic boundaries (harmonic potentials will keep solvents in the box)")
          if (nwords.ge.7) then
c           we can read the box dimensions from the option
            boxOrigin(1) = getDouble(words(2),ZERO)
            boxOrigin(2) = getDouble(words(3),ZERO)
            boxOrigin(3) = getDouble(words(4),ZERO)
            boxTop(1) = getDouble(words(5),ZERO)
            boxTop(2) = getDouble(words(6),ZERO)
            boxTop(3) = getDouble(words(7),ZERO)
            SolventBox = .true.
            SolventCap = .false.
          endif                                  
        else if (words(1).eq.'cap') then
          boundary = SPHERICALHARMONIC
          call printLine(INFO,"Requested the use of a spherical solvent cap")
          if (nwords.ge.6) then
c           we can read the cap coordinates/parameters fmro the option
            cap(1) = getDouble(words(2),ZERO)
            cap(2) = getDouble(words(3),ZERO)
            cap(3) = getDouble(words(4),ZERO)
            capradius = getDouble(words(5),ZERO)
            capconstant = getDouble(words(6),ZERO)
            SolventCap = .true.
            SolventBox = .false.
          endif
        else if (words(1).eq.'none') then
          boundary = NOBOUNDARIES
          call printLine(INFO,"Using no boundary conditions (gas-phase)")
        else if (words(1).eq.'solvent') then
          boundary = AGREESOLVENT
          call printLine(INFO,"Using whatever boundary conditions the solvent desires")
        else
          write(printstring,*) "Could not determine boundary conditions from ",words(1)
          call printLine(WARNING,printstring)
          call printLine(INFO,"Using whatever boundary conditions the solvent desires")
          boundary = AGREESOLVENT
        endif
      endif

c     get the non-bonded cutoff
      if (getOption('cutoff')) cut = getDouble(words(1),cut)
      if (cut.le.ZERO) cut = ZERO
      cut2 = cut**2
      write(printstring,*) "Non-bonded cutoff = ",cut," angstroms."
      call printLine(INFO,printstring)
      
c     get the feather factor
      if (getOption('feather')) feather = getDouble(words(1),feather)      
      if (feather.le.ZERO) feather = ZERO
      if (feather.gt.cut) feather = cut
      write(printstring,*) "Non-bonded cutoff feathered over the last ",feather," angstroms."
      call printLine(INFO,printstring)
      feather = (cut-feather)
      feather2 = feather**2
      featherfac = 1.0 / (cut2 - feather2)

c     now get the type of cutting
      if (getOption('cuttype')) then
        call lowercase(words(1))
        if (words(1).eq.'residue') then
          cuttype = RESIDUEMOLECULE
        else if (words(1).eq.'molecule') then
          cuttype = MOLECULEMOLECULE
        else
          write(printstring,*) "Could not interpret cuttype from ",words(1)
          call prettyLine(WARNING,printstring)
          call printLine(WARNING,"Using default cuttype instead")
        endif
      endif
      
      if (cuttype.eq.MOLECULEMOLECULE) then
        call printLine(INFO,"Using protein molecule/solute-molecule/solvent-molecule based cut-off.")
      else
        call printLine(INFO,"Using protein-residue/solute-molecule/solvent-molecule cut-off")
      endif
c
c     JM. Check if we want to run Surface Area Calculations
c
      if (getOption('surface')) then
        doSA = .true.
        probe = -1
        write(printstring,*) 'Doing Surface Area calculations'
        call printLine(INFO,printstring)        
        do j=1,nwords
          if (words(j)(1:1) .eq. '#') goto 91
          call lowercase(words(j))
          if (words(j) .eq. 'quality') then
            if (j .eq. nwords) then
              SAquality = 3
              write(printstring,*) 'SA quality is not determined ',
     &          ' assuming a value of ',SAquality
              call printLine(INFO,printstring)
            else
              SAquality = getInt(words(j+1),3)
              write(printstring,*) 'SA quality set to ',SAquality
              call printLine(INFO,printstring)
            endif
          else if (words(j) .eq. 'probe') then
            probe = getDouble(words(j+1),1.4d0)
            write(printstring,*) 'Probe radius set to ',probe
            call printLine(INFO,printstring)
          endif
        enddo
 91     continue
        if (SAquality .eq. -1) then
c         Problem, we did not read quality
          SAquality = 3
          write(printstring,*) 'SA quality is not determined assuming a value of ',SAquality
          call printLine(INFO,printstring)
        endif
c       Now specify the number of points on the surface of each sphere
        if (SAquality .eq. 1) then
          nsurfpoints = 12
        else if (SAquality .eq. 2) then
          nsurfpoints = 42
        else if (SAquality .eq. 3) then
          nsurfpoints = 162
        else if (SAquality .eq. 4) then 
          nsurfpoints = 642
        endif     
        if (probe .eq. -1) then
c         Problem, we dit not read a probe
          probe = 1.4d0
          write(printstring,*) 'Probe radius not determined assuming a value of ',probe
          call printLine(INFO,printstring)          
        endif
c       Set SurfCut
        SurfCut = 8.0d0
        SurfCut2 = SurfCut **2
      endif
c
c     JM Check if we want to run Generalised Born calculations
c
      if (getOption('born')) then
        doGB = .true.
        write(printstring,*) 'Doing Generalised Born calculations'
        call printLine(INFO,printstring)
        BornCut = -1.0d0
        gbproteins = .false.
        gbthreshold = -1.0d0
        do j=1,nwords
          if (words(j)(1:1) .eq. '#') goto 92
          call lowercase(words(j))
          if (words(j) .eq. 'cut') then
            BornCut = getDouble(words(j+1),25.0d0)
            write(printstring,*) 'BornCut set to ',BornCut
            call printLine(INFO,printstring)             
          endif
          if (words(j) .eq. 'threshold') then
            gbthreshold = getDouble(words(j+1),0.01d0)
            write(printstring,*) 'GBthreshold set to ',gbthreshold
            call printLine(INFO,printstring)           
          endif
          if (words(j) .eq. 'proteins') then
            gbproteins = .true.
            write(printstring,*) 'Using Born corrections for Proteins'
            call printLine(INFO,printstring)  
          endif
        enddo
 92     continue
        if (BornCut .lt. ZERO) then
          BornCut = 25.0d0
          write(printstring,*) 'Default BornCut set to ',BornCut
          call printLine(INFO,printstring)             
        endif
        BornCut2 = BornCut**2
      endif

      if (getOption('softcoreprotein')) then
        ProteinSoftCore = getLogical(words(1), ProteinSoftCore)

        if (ProteinSoftCore) then
          call prettyLine(INFO, "Softening residue-residue interactions")
        endif
      endif

c
c     JM
c
      if (getOption('softcoreparams')) then
        do j=1,nwords
          if (words(j) .eq. 'coul') CoulSoftPower = getInt(words(j+1),1)       
          if (words(j) .eq. 'lj') LJSoftPower = getInt(words(j+1),1)       
          if (words(j) .eq. 'gb') GBSoftPower = getInt(words(j+1),1)       
          if (words(j) .eq. 'delta') DeltaCore = getDouble(words(j+1),0.5d0)
          if (words(j) .eq. 'deltacoul') DeltaCoreCoulL = getDouble(words(j+1),2.0d0)
c       Have not yet implemented general soft-core
c          if (words(j) .eq. 'power') SoftPower = getDouble(words(j+1),6.0d0)
c         Different types of soft-core potentials
          if (words(j) .eq. 'amber') SoftType = SOFTAMB
          if (words(j) .eq. 'old') SoftType = SOFTOLD
          if (words(j) .eq. 'soft66') SoftType = SOFT66
        enddo
        if (SoftType.eq.SOFT66) then
          SoftPower = 6.0d0
          DeltaCoreCoul = DeltaCoreCoulL**6
        else if (SoftType.eq.SOFTAMB) then
          DeltaCoreCoul = DeltaCoreCoulL
        endif
        SoftPowerInv = ONE / SoftPower
c       SG these two are not used since general soft-core is not yet supported
        SoftPowerInv6 = 6.0d0 / SoftPower
        SoftPower2 = SoftPower / 2.0d0
        write(printstring,*) 'SoftCore Parameters for this simulation :'
        call printLine(INFO,printstring)          
        write(printstring,*) 'Coul ',CoulSoftPower,' LJ ',LJSoftPower,
     .    ' GB ',GBSoftPower,' Delta ',DeltaCore, ' DeltaCoul ',DeltaCoreCoulL, ' ( ',
     .     DeltaCoreCoul,' ) Power ',int(SoftPower)
        call printLine(INFO,printstring)
        if (SoftType.eq.SOFTAMB) then
          call printLine(INFO,"Using Amber-like (2-6) soft-core potential")
        else if (SoftType.eq.SOFT66) then
          call printLine(INFO,"Using 6-6 soft-core potential")       
        else
          call printLine(INFO,"Using old-type soft-core potential")       
        endif
      endif
c
c     JM 
c
c$$$      if (getOption('quickpotential')) then
c$$$        do j=1,nwords
c$$$          if (words(j) .eq. 'gboff') OnGBQuickPot = .false.
c$$$          if (words(j) .eq. 'gbon') OnGBQuickPot = .true.
c$$$          if (words(j) .eq. 'dddon') OnDDDQuickPot = .true.
c$$$          if (words(j) .eq. 'dddoff') OnDDDQuickPot = .false.
c$$$        enddo
c$$$        write(printstring,*) 'QuickPotential Parameters :'
c$$$        call printLine(INFO,printstring)
c$$$        write(printstring,*) ' GB On ? ',OnGBQuickPot,' DDD On ? ',OnDDDQuickPot
c$$$        call printLine(INFO,printstring)
c$$$      endif

c     get the names of all of the parameter files and read them in
c     (as the parameter files may change some of the default options)
c     Keep reading parameter files while we have them
      do i=1,nkeys
c       7 in the line below stands for the lenght of 'dump'
        if (keys(i)(1:7).eq.'parfile') then
          call split(values(i))
c         copy words(1) to filename as readParFile will overwrite 'words' when
c         it splits any lines    
          filename = words(1)
          call readParFile(filename)
        endif
      enddo

c     get the names of all of the protein files and read them in
      i = 0
101   continue
        i = i+1
        call writeNumberString('protein',i,tempstring)
        if (getOption(tempstring)) then
          filename = words(1)
          call loadProtein(filename)
          goto 101
        endif
c     end of do-while

c     SG
c     Try to read in number of solute theta kinds
      NThetaSolKind = 0
      if (getOption('nthetasolutes').or.getOption('nfrag')) then
        NThetaSolKind = getInt(words(1),0)  
        write(printstring,*)"Will consider the first ",NThetaSolKind, "soluteX as theta-solutes"
        call printLine(INFO,printstring)
      endif      
      nkind = 0

c     get the names of all of the solute files and read them in
      i = 0
102   continue      
        i = i+1
        call writeNumberString('solute',i,tempstring)
        if (getOption(tempstring)) then
          filename = words(1)
          nold = NSolutes
          call loadSolute(filename)
c         SG, check if we should tag this solute as a theta solute
          if (i.le.NThetaSolKind) then
            nkind = nkind + 1
c           Set default values of fragment variables            
            do j=nold+1,NSolutes
              ThetaKind(j) = nkind
              SoluteTheta(j) = ZERO
              SoluteHydScale(j) = ONE
            enddo
            NThetaSolutes = NThetaSolutes + NSolutes-nold
          else
            do j=nold+1,NSolutes
              ThetaKind(j) = 0
              SoluteTheta(j) = ONE
              SoluteHydScale(j) = ONE
            enddo
          endif
          goto 102
        endif        
c     end of do-while

c     SG, set initial theta values
      if (getOption('theta0')) then
        do i=1,min(NThetaSolutes,nwords)
          SoluteTheta(i) = getDouble(words(i),SoluteTheta(i)) 
        enddo
      endif 

c     SG, set hydration scaling variable
      if (getOption('hydscaling')) then
        do i=1,min(NThetaSolutes,nwords)
          SoluteHydScale(i) = getDouble(words(i),SoluteHydScale(i)) 
        enddo
      endif 

      UseConcCorr = .true.
      if (getOption('conccorr')) then
        UseConcCorr = (words(1).eq."on") 
      endif
      if (UseConcCorr) then
        write(printstring,*)"Using concentration correction"
        call printLine(INFO,printstring)
      else
        write(printstring,*)"NOT using concentration correction"
        call printLine(INFO,printstring)
      endif

c     get the names of all of the gc- solute files and read them in
      i = 0
103   continue
        i = i+1
        call writeNumberString('grand',i,tempstring)
        if (getOption(tempstring)) then
          filename = words(1)
          call loadgc(filename)
          goto 103
        endif
c     end of do-while

c     get the names of all of the solvent files and read them in
      i = 0
104   continue
        i= i+1
        call writeNumberString('solvent',i,tempstring)
        if (getOption(tempstring)) then
          filename = words(1)
          call loadSolvent(filename)
          goto 104
        endif
c     end of do-while

      write(printstring,*) "Simulation consists of ",NProteins," protein chains, ",NSolutes,
     .                     " solute molecules and ",NSolvents," solvent molecules"
      call printLine(INFO,printstring)

c     check that at least something will be loaded!
      if (NProteins.eq.0 .and. NSolutes.eq.0 .and. NSolvents.eq.0) then
        call closeProgram(-1,"Nothing has been loaded - ProtoMC cannot continue!")
      endif

c     Assign the Surface Area terms to each CLJ type that was loaded
      if (doSA) call assignSurface
c     Assign Generalised Born terms to each CLJ type that was loaded
      if (doGB) call assignBorn

c     now assign templates etc. for all of the loaded molecules
      call assignSystem

      if (NSolutes.gt.0) then
        do i=1,NSolutes
c         trying to get possible solute groups
          call writeNumberString('group',i,tempstring)
          if (getOption(tempstring)) then
c           get the number of solutes in this group
            NGroupSols(i) = getInt(words(1),0)
c           loop over the amount of solutes in a group
            do j=1,getInt(words(1),0)
              SolGroup(i,j) = getInt(words(j+1),0)
            enddo
c           Here I'm setting the number of groups to the number of groups user-defined.
c           I'm not sure whether it should be like this, or it should include the automatic.
c           But in that case, it is simply equal to the number of solutes...
            NSolGroups = i
          else
            NGroupSols(i) = 1
            SolGroup(i,1) = i
          endif
        enddo
      endif
      
      if (NSolutes.gt.0) then
c       get the number of dual topology solutes - do this here as we may need to 
c       modify the dummies and this can only be done after the solutes have been loaded
        i = 0
200     continue
          i = i+1
          call writeNumberString('dualtopology',i,tempstring)
          if (getOption(tempstring)) then
c           format of the line is 'dualtopology lowgroup highgroup (syncdummy) (synctrans) (syncrot) #comment'
            if (nwords.lt.2) then
              call printLine(WARNING,"dualtopology line needs two parameters (lowgroup and highgroup)")
              i = i-1
              goto 200
            endif
            
c           read in both solutes
            lowgroup = getInt(words(1),0)
            highgroup = getInt(words(2),0)                        
          
            if (lowgroup.gt.0 .and. lowgroup.le.NSolutes .and.
     .          highgroup.gt.0 .and. highgroup.le.NSolutes .and. 
     .          lowgroup.ne.highgroup) then
          
              lgth0 = lineLength(SoluteName(lowgroup))
              lgth1 = lineLength(SoluteName(highgroup))
              
              write(printstring,*) "Dual topology simulation with group ",
     .              lowgroup," at lambda=0, and group ",
     .              highgroup," at lambda=1. If groups are not defined, ",
     .              "each group corresponds to the solute with the same number."
              call printLine(INFO,printstring)
            
              if (DualTopologySolute(SolGroup(lowgroup,1)).ne.0.and.
     .        DualTopologySolute(SolGroup(lowgroup,NGroupSols(lowgroup))).ne.0) then
                write(printstring,*) "Group ",
     .               lowgroup," already has a dual topology partner!"
                call printLine(WARNING,printstring)
                i = i-1
                goto 200
              else if (DualTopologySolute(SolGroup(highgroup,1)).ne.0.and.
     .        DualTopologySolute(SolGroup(highgroup,NGroupSols(highgroup))).ne.0) then
                write(printstring,*) "Group ",
     .               highgroup," already has a dual topology partner!"
                call printLine(WARNING,printstring)
                i = i-1
                goto 200
              endif
            
              do j=1,NGroupSols(lowgroup)
                DualTopologySolute(SolGroup(lowgroup,j)) = highgroup
                DualTopologySyncTrans(SolGroup(lowgroup,j)) = .false.
                DualTopologySyncRot(SolGroup(lowgroup,j)) = .false.
                DualTopologySyncInt(SolGroup(lowgroup,j)) = .false.
              enddo
              do j=1,NGroupSols(highgroup)
                DualTopologySolute(SolGroup(highgroup,j)) = -lowgroup
                DualTopologySyncTrans(SolGroup(highgroup,j)) = .false.
                DualTopologySyncRot(SolGroup(highgroup,j)) = .false.
                DualTopologySyncInt(SolGroup(highgroup,j)) = .false.
              enddo
              
c             loop through the rest of the string to see if we sync moving or 
c             sync dummies
              if (nwords.eq.2) goto 200
              
              do j=3,nwords
                if (words(j)(1:1).eq.'#') goto 201
                call lowercase(words(j))
                if (words(j).eq.'syncdummy') then
c                 we will sync the solute dummies - this ensures that the two
c                 solutes are rotated about the same point, thus ensuring that
c                 the solutes stay together
                  call printLine(INFO,"Dummy atoms of both solutes will be synched")
                  call syncSoluteDummies(lowgroup,highgroup)
                else if (words(j).eq.'synctrans') then
c                 we will translate these solutes together
                  do k=1,NGroupSols(highgroup)
                    DualTopologySyncTrans(SolGroup(highgroup,k)) = .true.
                  enddo
                  do k=1,NGroupSols(lowgroup)
                    DualTopologySyncTrans(SolGroup(lowgroup,k)) = .true.
                  enddo
                  call printLine(INFO,"Solutes will be translated together")
                else if (words(j).eq.'syncrot') then
c                 we will rotate these solutes together
                  do k=1,NGroupSols(highgroup)
                    DualTopologySyncRot(SolGroup(highgroup,k)) = .true.
                  enddo
                  do k=1,NGroupSols(lowgroup)
                    DualTopologySyncRot(SolGroup(lowgroup,k)) = .true.
                  enddo
                  call printLine(INFO,"Solutes will rotated together")
                endif                  
              enddo
201           continue
            else
              write(printstring,*) "Problem specifying dual topology solutes, ",
     .           "lowgroup = ",lowgroup," and highgroup = ",highgroup
              call printLine(WARNING,printstring)
              i = i-1
              goto 200
            endif
            goto 200
          endif
c       end of do-while loop
      endif  

c     Check wheter or not a soft core is to be used
c     The format is softcore1 solute 1, softcore2 protein 2, softcore3 solvent 45
c     A softcore is applied if either atom is part of a softcore residue
      i = 0
 210  continue
      i = i+1
      call writeNumberString('softcore',i,tempstring)
      if (getOption(tempstring)) then 
        if (nwords .lt. 2) then
          call printLine(WARNING,"softcore format is <moltype> <molnum>")
          i = i-1
          goto 210
        endif
        call lowercase(words(2))
        if (words(2) .eq. 'all') then 
          imol = -1
        else
          imol = getInt(words(2),0)
        endif
        if (imol .eq. 0) then
          call printLine(WARNING,"cannot get molecule number !")
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)            
        endif
        call lowercase(words(1))
        if (words(1) .eq. 'solute') then
          if (imol .eq. -1) then
            do j=1,NSolutes
              SoluteSoftCore(j) = .true.
            enddo
          else
            SoluteSoftCore(imol) = .true.
          endif
          write(printstring,*) 'Turning softcore on for Solute ',words(2)  
        else if (words(1) .eq. 'gcsolute') then
         if (imol .eq. -1) then
            do j=1,NSolutes
              GCsoluteSoftCore(j) = .true.
            enddo
          else
            GCsoluteSoftCore(imol) = .true.
          endif
          write(printstring,*) 'Turning softcore on for GCsolute ',words(2)                
        else if (words(1) .eq. 'protein') then
          write(printstring,*) 'No Softcore for Proteins ',words(2)
        else if (words(1) .eq. 'solvent') then
          write(printstring,*) 'No Softcore for Solvent ',words(2)
        else
          call printLine(WARNING,"cannot interpret type of molecule !")
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)            
          i = i-1
          goto 210
        endif

        call printLine(INFO,printstring)
        goto 210
      endif
c     end of do-while loop
      

c     get the temperature
      if ((.not.doTemperatureRE).and.getOption('temperature')) then
        temperatureC = getDouble(words(1),temperatureC)
        temperature = temperatureC + DEG2KEL
      elseif (doTemperatureRE) then
        temperatureC = temval(mytemperature)
        write(6,*)myid," My Just Stored temperature ",temperatureC
        temperature = temperatureC + DEG2KEL
        if (doLambdaRE) then
          pt_beta_t(mylambda,mytemperature) = dble(1.0) / (GASR*(temperature))
        endif
      endif
      write(printstring,*) "Simulation temperature = ",temperatureC," C, (",
     .                              temperature," K)"
      call printLine(INFO,printstring)
      mt_beta_t = dble(1.0) / (GASR*temperature)
      write(6,*) myid,"GetOptions temperaturec__",temperatureC

c     get the REST2 effective temperature
      resttempC = temperatureC
      if (getOption('resttemp')) resttempC = getDouble(words(1),temperatureC)
      resttemp = resttempC + DEG2KEL
      rest_beta = dble(1.0) / (GASR*resttemp)
      restscale1 = rest_beta / mt_beta_t
      restscale2 = (mt_beta_t+rest_beta)/(2.0d0*mt_beta_t)
      restscale3 = 1.0d0 

c     get the B factor for GCMC simulations
      if ((.not.doMultiGCMC).and.getOption('potential')) then
         B = getDouble(words(1),B)
        write(printstring,*) "Simulation B value / Adams parameter = ",B," AU"
        call printLine(INFO,printstring)
      endif
      if (doMultiGCMC) then
        B = adamparam(myadam)
        write(printstring,*) "Simulation B value / Adams parameter = ",B," AU"
        call printLine(INFO,printstring) 
      endif

c     Get the JAWS biasing term for JAWS-II simulations
      if ((.not.doMultiJaws2).and.getOption('jbias')) then
         bias = getDouble(words(1),bias)
        write(printstring,*) "JAWS biasing term = ",bias," kcal/mol"
        call printLine(INFO,printstring)
      endif
      if (doMultiJaws2) then
        bias = alljbias(myjbias)
        write(printstring,*) "JAWS biasing term = ",bias," kcal/mol"
        call printLine(INFO,printstring) 
      endif

      if (getOption('x').or.getOption('lenx')) then
      gridx = getDouble(words(1),gridx)
      write(printstring,*) "Grid size in X dimension" ,gridx
      call printLine(INFO,printstring)
      endif
      if (getOption('y').or.getOption('leny')) then
      gridy = getDouble(words(1),gridy)
      write(printstring,*) "Grid size in Y dimension" ,gridy
      call printLine(INFO,printstring)
      endif
      if (getOption('z').or.getOption('lenz')) then
      gridz = getDouble(words(1),gridz)
      write(printstring,*) "Grid size in Z dimension" ,gridz
      call printLine(INFO,printstring)
      endif

      if (getOption('spacing')) then
      space = getDouble(words(1),space)
      write(printstring,*) "Grid granularity" ,space, "A"
      call printLine(INFO,printstring)
      endif

c     get grid centers
      if (getOption('centerx')) then
      centerx = getDouble(words(1),centerx)
      write(printstring,*) "Grid center in X dimension" ,centerx
      call printLine(INFO,printstring)
      coordx = centerx - gridx / 2.0
      write(printstring,*) "Grid origin in X dimension" ,coordx
      call printLine(INFO,printstring)
      endif

      if (getOption('centery')) then
      centery = getDouble(words(1),centery)
      write(printstring,*) "Grid center in Y dimension" ,centery
      call printLine(INFO,printstring)
      coordy = centery - gridy / 2.0
      write(printstring,*) "Grid origin in Y dimension" ,coordy
      call printLine(INFO,printstring)
      endif

      if (getOption('centerz')) then
      centerz = getDouble(words(1),centerz)
      write(printstring,*) "Grid center in Z dimension" ,centerz
      call printLine(INFO,printstring)
      coordz = centerz - gridz / 2.0
      write(printstring,*) "Grid origin in Z dimension" ,coordz
      call printLine(INFO,printstring)
      endif

c     get grid origins
      if (getOption('originx')) then
       coordx = getDouble(words(1),coordx)
      write(printstring,*) "Grid origin in X dimension" ,coordx
      call printLine(INFO,printstring)
      endif

      if (getOption('originy')) then
       coordy = getDouble(words(1),coordy)
      write(printstring,*) "Grid origin in Y dimension" ,coordy
      call printLine(INFO,printstring)
      endif
c"
      if (getOption('originz')) then
       coordz = getDouble(words(1),coordz)
      write(printstring,*) "Grid origin in Z dimension" ,coordz
      call printLine(INFO,printstring)
      endif

      if (getOption('gcmc')) then
      GCtheta = getDouble(words(1),GCtheta)
      write(printstring,*) "Running a GCMC sim",GCtheta
      call printLine(INFO,printstring)
        do i=1,NGCSolutes
         thetaparam(i) = GCtheta
        enddo
      endif

      if (getOption('jaws2')) then
      jaws2 = getDouble(words(1),jaws2)
      write(printstring,*) "Running a JAWS-2 sim"
      call printLine(INFO,printstring)
        do i=1,NGCSolutes
        thetaparam(i) = jaws2
        enddo
      endif

      if (getOption('thres')) then
      thres = getDouble(words(1),thres)
      write(printstring,*) "JAWS threshold for theta sampling is" ,thres
      call printLine(INFO,printstring)
      endif


      if (getOption('jaws1')) then
      jaws1 = getDouble(words(1),jaws1)
      write(printstring,*) "Running a JAWS-1 sim"
      call printLine(INFO,printstring)
        do i=1,NGCSolutes
        thetaparam(i) = jaws1
        enddo
      endif

c     read perturbeIntra flag
      perturbeIntra = .true.
      if (getOption('perturbeintra')) then
        call lowercase(words(1))
        if (words(1).eq.'off') then
          perturbeIntra = .false.
          write(printstring,*) "Will not perturbe intramolecule solute energies"
          call printLine(INFO,printstring)
        endif
      endif 

c     see whether we are doing an npt simulation, then get the pressure if we are
      if (NSolvents.gt.0) then
        if (getOption('pressure')) then
          pressure = getDouble(words(1),ZERO)
          if (pressure .ne. ZERO) then

            intPressure = pressure*P2INTERNALP
c           get the maximum change in volume
            maxVolChange = dble(NSolvents)*0.1D+00
            if (getOption('maxvolchange')) 
     .           maxVolChange = getDouble(words(1),maxVolChange)

c           must have periodic boundaries or spherical harmonic
            if (boundary .eq. PERIODICBOUNDARIES) then
              write(printstring,*) "Running an NPT periodic boundaries ",
     .           "simulation with pressure = ",
     .           pressure," atm (",pressure*ATM2PASCAL," pascals) and with ",
     .           "a maximum volume change of ",maxVolChange," A^3"
              call prettyLine(INFO,printstring)
           else if (boundary .eq. SPHERICALHARMONIC) then
c            calculate a rough 'volume' of this cap. This will be the volume
c            of the cap, plus an additional radius of kT into the harmonic potential
             volume = sqrt(GASR*temperature / capconstant)

             write(printstring,*) "Solvent cap volume calculated using extra radius of ",
     .                            volume, " A "
             call printLine(INFO,printstring)

             volume = (4.0/3.0) * PI * (capradius+volume)**3

             write(printstring,*) "Solvent cap volume = ",volume," A^3"
             call printLine(INFO,printstring)

             write(printstring,*) "Running an NPT solvent cap simulation ",
     .           "with pressure = ",pressure," atm (",pressure*ATM2PASCAL,
     .           " pascals) and with a maximum volume change of ",
     .           maxVolChange," A^3"
             call prettyLine(INFO,printstring)

           write(printstring,*)"Adjusting the volume when using a cap is not recommended. This may seriously harm your simulation!"
             call printLINE(WARNING,printstring)
            else
              call printLine(WARNING,"Must have periodic boundaries for an NPT simulation!")
              call prettyLine(INFO,"Running an NVT simulation")
              pressure = ZERO
              intPressure = ZERO
              maxVolChange = ZERO
            endif
          endif
        else
          call prettyLine(INFO,"Running an NVT simulation")
          pressure = ZERO
          intPressure = ZERO
          maxVolChange = ZERO
        endif
      endif
                         
c-----------------------------------------------------
c     ####### SIMULATION CONTROL PARAMETERS ########
c     ####### Some of these may already have #######
c     ####### been set by the parfiles  ############
c-----------------------------------------------------
      
c     Now whether or not we are doing preferential sampling,
c     and the value of the pref sampling parameter
      if (getOption('prefsampling')) then
        svnPrefSol = getInt(words(1),0)
        if (svnPrefSol.le.0 .or. svnPrefSol.gt.NSolutes) then
          call printLine(INFO,"Turning off preferential sampling")
          prefSample = .false.
        else
          prefSample = .true.
c         see if the line also contains the pref sampling parameter          
          if (nwords.ge.2) svnProbParam = getDouble(words(2),svnProbParam)
        
          lgth = lineLength(SoluteName(svnPrefSol))
          write(printstring,*) "Performing preferential sampling with ",
     .         "sampling parameter = ",svnProbParam," and centered ",
     .         "on solute ",SoluteName(svnPrefSol)(1:lgth)," (",
     .         svnPrefSol,")"
          call printLine(INFO,printstring)
        endif
      endif

      if (NSolvents.gt.0 .and. .not.prefSample) 
     .          call printLine(INFO,"Not performing preferential sampling")
      

c##########################      
c     finally get the random number seed - the default value is based on the current time
      call getDateAndTime(year,month,day,hour,mins,secs)
      
c     generate a different seed based on this time      
      ranseed = 1+hour+secs*100+mins*10000+day*1000000+month*100000000
      
      if (getOption('ranseed')) ranseed = getInt(words(1),ranseed)
      
      write(printstring,*) "Got initial random number seed ",ranseed
      call printLine(INFO,printstring)
      
c     make into an i6 or i7 positive, odd integer...
      if (ranseed.lt.0) ranseed = -ranseed
      if (mod(ranseed,2).eq.0) ranseed = ranseed - 1
      if (ranseed.lt.100000) ranseed = ranseed+110000
      write(printstring,*) "Sanitised random number seed = ",ranseed
      call prettyLine(INFO,printstring)
      
c     initialise the random number generator with this seed
      call sgrnd(ranseed)      

c     since we have the time, why not print out the time that the simulation will start!
      write(printstring,10) hour,mins,secs,day,month,year
10    format("Starting simulation at ",i2,":",i2.2,":",i2.2," on the ",i2,"/",i2,"/",i4)      
      call printLine(INFO,printstring)
      
      end
