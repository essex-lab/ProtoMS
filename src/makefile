#-----------------------------------------------
#  Makefile for ProtoMS
#     Christopher Woods
#       October 2004
#-----------------------------------------------

# Do you want to compile the RETI handler? (HAVE YOU GOT MPI?)
HAVE_MPI = 1

# What is the name of the fortran compiler that you want to use?
fortrancompiler = gfortran
mpi_fortrancompiler = mpif77

# The GNU gfortran compiler
#fortrancompiler = gfortran
#fortrancompiler = g77
#fortrancompiler = pgf77
#fortrancompiler = ifort

# Is openmp available? It needs to be to use the reti openmp driver
HAVE_OPENMP = 1

# The flag used by the fortran compiler to turn on openmp
#openmp_flags = -fopenmp

# What compilation flags do you want to use? Note the you need to turn
# on 'long lines'. For g77 you do this via the '-ffixed-line-length-none' option,
# while the Portland Group compilers use the '-Mextend' option.
# Either modify the line below, or uncomment the line that matches your
# compiler / architecture
fflags = -m64 -O3 -I. -fomit-frame-pointer -ffixed-line-length-none

#Compile flags for g77
#fflags = -O2 -malign-double -ffixed-line-length-none -Wall
#Compile flags for g77 with gprof profiling and debugging
#fflags = -g -pg -O2 -malign-double -ffixed-line-length-none -Wall -ffortran-bounds-check
#Compile flags for g95
#fflags = -O2 -I. -malign-double -Wall -ffixed-line-length-132
#Compile flags for pgf77
#fflags = -O2 -fast -Mextend
#Compile flags for ifort
#fflags = -O3 -132

# the link flags - either edit the below line or uncomment the line that matches
# your system
#LFLAGS = -static
# link flags for static executable
#LFLAGS = -static
# link flags for dynamic executable
#LFLAGS = 

# If you do not have python2 installed then uncomment the below line
# (you will know that you don't have python2 installed if the compilation
# fails with the error 'getid.py: Command not found', or something similar
#NOPYTHON2 = 1

#*******************************************************
#  This is the list of all source files used to
#  build ProtoMS. If you want to add a new file
#  to ProtoMS then just add it to this list.
#*******************************************************

OBJS = addacid.o \
       addangle.o \
       addbase.o \
       addbond.o \
       addbondharmonicrestraint.o \
       addcartesianharmonicrestraint.o \
       addcartesianflatbottomrestraint.o \
       adddihedralharmonicrestraint.o \
       addchainangle.o \
       addchainbond.o \
       addchaindihedral.o \
       addchainureybradley.o \
       addclj.o \
       adddihedral.o \
       addgcsoluteaxes.o \
       addid.o \
       addmolangle.o \
       addmolatom.o \
       addmolbond.o \
       addmoldihedral.o \
       addmolureybradley.o \
       addproteinangle.o \
       addproteinbond.o \
       addproteindihedral.o \
       addproteinureybradley.o \
       addsoluteangle.o \
       addsolutebond.o \
       addsoluteclj.o	\
       addsolutedihedral.o \
       addsolutedummies.o \
       addsolutenb.o \
       addsoluteureybradley.o \
       addsolventaxes.o \
       addureybradley.o \
       amber.o \
       angleenergy.o \
       assignborn.o \
       assigngcsolutetemplate.o \
       assigngcsolutetemplates.o \
       assignresidueoverrides.o \
       assignresiduetemplate.o \
       assignrotamerlibrary.o \
       assignsolutetemplates.o \
       assignsolutevariables.o \
       assignsolventtemplate.o \
       assignsolventtemplates.o \
       assignsurface.o \
       assignsystem.o \
       atomname.o \
       biaspotential.o \
       bondedresresenergy.o \
       bondenergy.o \
       bndharmonicrestraintenergy.o \
       buildprotein.o \
       buildresidue.o \
       buildsolute.o \
       buildsolventbox.o \
       calcangle.o \
       calcaverageenergy.o \
       calcaveragegbenergy.o \
       calcaveragesaenergy.o \
       calcaverageextraenergies.o \
       calcbond.o \
       calcdihedral.o \
       calcgeneraldihedral.o \
       calculatecljpairs.o \
       calculatesolventbox.o \
       calculatetotalsurface.o \
       calcureybradley.o \
       changelambda.o \
       charge.o \
       checkbonded.o \
       checkgottotal.o \
       cleansoluteparameters.o \
       cleargcsolutecoordstack.o \
       clearid.o \
       clearmovelist.o \
       clearresiduecoordstack.o \
       clearresidueenergystack.o \
       clearsavesolutegb.o \
       clearsavesolventgb.o \
       clearsaveresiduegb.o \
       clearsolutecoordstack.o \
       clearsolventcoordstack.o \
       clearsplitmovelist.o \
       clearsurfaces.o \
       cljid.o \
       closefile.o \
       closefiles.o \
       closeprogram.o \
       convertres.o \
       copysolvent.o \
       covalentradius.o \
       createconnectivitymap.o \
       createresresconnectivitymap.o \
       cross.o \
       deg2rad.o \
       deletionmove.o \
       deprotonate.o \
       detailedresidueenergy.o \
       detailedsoluteenergy.o \
       dgenergies.o \
       dihedralenergy.o \
       dihedralharmonicrestraintenergy.o \
       dihenergy.o \
       dot.o \
       dualtopologysasa.o \
       dump.o	\
       duplicatesolute.o \
       elementname.o \
       epslon.o \
       equilibrate.o \
       extraenergy.o \
       fakesimulate.o \
       findangleindex.o \
       findangleparam.o \
       findatom.o \
       findbondindex.o \
       findbondparam.o \
       findchainangle.o \
       findchainatom.o \
       findchainbond.o \
       findchaindihedral.o \
       findchaintemplate.o \
       findchainureybradley.o \
       finddihedralindex.o \
       finddihedralparam.o \
       findgb.o \
       findgcsolutetemplate.o \
       findgcsolutetemplateatom.o \
       findmoleculetemplate.o \
       findmoltemplateangle.o \
       findmoltemplateatom.o \
       findmoltemplatebond.o \
       findmoltemplatedihedral.o \
       findmoltemplateureybradley.o \
       findoverride.o \
       findresidue.o \
       findresiduetemplate.o \
       findSA.o \
       findsolresidue.o \
       findsoluteatom.o \
       findsolventatom.o \
       findsolventtemplate.o \
       findsolventtemplateatom.o \
       findtemplate.o \
       findtemplateangle.o \
       findtemplateatom.o \
       findtemplatebond.o \
       findtemplatedihedral.o \
       findtemplateureybradley.o \
       findureybradleyindex.o \
       findureybradleyparam.o \
       fixresidues.o	\
       flagmovedbits.o \
       flatbottomrestraintenergy.o \
       fliphistidine.o \
       fullsplit.o \
       gbenergy.o \
       gbsavesystem.o \
       gcsolutebulkenergy.o \
       gcsoluteenergy.o \
       gcsolutegcsoluteenergy.o \
       gcsolutegcsolutesingleenergy.o \
       gcsolutemove.o \
       gcsolutesinglesolventenergy.o \
       gcsolutesolventenergy.o \
       genproteinatom.o \
       gensoluteatom.o \
       getangle.o \
       getatomname.o \
       getaxisvectors.o \
       getbond.o \
       getchaintemplate.o \
       getclj.o \
       getcoordinates.o \
       getdateandtime.o \
       getdihedral.o \
       getdihwells.o \
       getdofinfo.o \
       getdouble.o \
       geteffectivebornradii.o \
       getgbiborn.o \
       getgbscale.o \
       gethighscale.o \
       getint.o \
       getlastresidue.o \
       getlogical.o \
       getlowscale.o \
       getmajorminoraxis.o \
       getmoveinfo.o \
       getmovetype.o \
       getnonbonded.o \
       getoption.o \
       getoptions.o \
       getperpendicular.o \
       getrange.o \
       getresiduecentreradii.o \
       getresiduecentreradius.o \
       getresintranonbonded.o \
       getresresbondindex.o \
       getresresinternonbonded.o \
       getsmallestdelta.o \
       getsolatomname.o \
       getsolutecentreradius.o \
       getsolutenonbonded.o \
       getsolutezmatrix.o \
       getureybradley.o \
       getvectoronsphere.o \
       getzmatdihedral.o \
       grnd.o \
       hardwallenergy.o \
       harmonicrestraintenergy.o \
       hasangle.o \
       hasatommoved.o \
       hasbond.o \
       hasdihedral.o \
       hasresiduemoved.o \
       hassoluteangle.o \
       hassolutebond.o \
       hassolutedihedral.o \
       hassolutemoved.o \
       hassoluteureybradley.o \
       hassolventmoved.o \
       hasureybradley.o \
       initialise.o \
       initialisearrays.o \
       initialiseextraenergies.o \
       initialiseprintsystem.o \
       initialiseprotein.o \
       initialisesolute.o \
       inresiduesplitmoved.o \
       insertionmove.o \
       insolutesplitmoved.o \
       insolventsplitmoved.o \
       interpretposition.o \
       interresidueresidueenergy.o \
       intraresidueproteinenergy.o \
       isbackbone.o \
       join.o \
       join2.o \
       lambdare.o	\
       lreenergy.o	\
       lambdamove.o \
       leftjustify.o \
       linelength.o \
       loadgc.o \
       loadprotein.o \
       loadsolute.o \
       loadsolvent.o \
       lowercase.o \
       main.o \
       makevariableangle.o \
       makevariablebond.o \
       makevariabledihedral.o \
       makevariablesolangle.o \
       makevariablesolbond.o \
       makevariablesoldihedral.o \
       markgbupdate.o \
       marksaupdate.o \
       metropolis.o \
       metropolis2.o \
       metropolis3.o \
       minsolventgcsolutedist2.o \
       minsolventproteindist2.o \
       minsolventsolutedist2.o \
       moveacid.o \
       movegcsolute.o \
       moverotamer.o \
       movesidechain.o \
       movesolute.o \
       movesolvent.o \
       movevolume.o \
       multiplemovesenergy.o \
       namestring.o \
       notbonded.o \
       openfile.o \
       openstreamfiles.o \
       parseargs.o \
       parsefile.o \
       perturbforcefield.o \
       popallresidueenergy.o \
       popdeletesolventcoords.o \
       popgcsolutecoords.o \
       poppdacont.o \
       popresiduecoords.o \
       popsolutecoords.o \
       popsolventcoords.o \
       positionstring.o \
       prettyline.o \
       printattemptaccept.o \
       printaverages.o \
       printaveragesgbenergy.o \
       printaveragesaenergy.o \
       printaverageextraenergies.o \
       printboxbar.o \
       printboxtext.o \
       printenergies.o \
       printline.o \
       printparameters.o \
       printpdb.o \
       printpdbatom.o \
       printpdbgcsolute.o \
       printpdbprotein.o \
       printpdbresidue.o \
       printpdbsolute.o \
       printpdbsolvent.o \
       processgeometryline.o \
       processidline.o \
       processhardwallline.o \
       processrestraintline.o \
       processsolventline.o \
       proteingcsenergy.o \
       proteinintraenergy.o \
       proteinproteinenergy.o \
       proteinresresbondlist.o \
       proteinsoluteenergy.o \
       proteinsolventenergy.o \
       protonate.o \
       protonnum.o \
       pushbackgcsolutecoords.o \
       pushbacknewresidueenergy.o \
       pushbackoldresidueenergy.o \
       pushbackpdacont.o \
       pushbackresiduecoords.o \
       pushbacksolutecoords.o \
       pushbacksolventcoords.o \
       rad2deg.o \
       ranint.o \
       ranrange.o \
       readang.o \
       readbnd.o \
       readborn.o \
       readchaintemplate.o \
       readchaintemplateangleline.o \
       readchaintemplateatomline.o \
       readchaintemplatebbatomline.o \
       readchaintemplatebondline.o \
       readchaintemplatedihedralline.o \
       readchaintemplateubline.o \
       readchaintemplatezmatline.o \
       readclj.o \
       readdih.o \
       readdouble.o \
       readgcsolutetemplate.o \
       readgcstemplateatomline.o \
       readgcstemplateinfoline.o \
       readinfo.o \
       readinfoline.o \
       readint.o \
       readlogical.o \
       readoverride.o \
       readparfile.o \
       readpdbline.o \
       readresidueoverride.o \
       readresiduetemplate.o \
       readrestart.o \
       readrotamer.o \
       readsurface.o \
       readrestemplateangleline.o \
       readrestemplateatomline.o \
       readrestemplatebackboneline.o \
       readrestemplatebondline.o \
       readrestemplatedihedralline.o \
       readrestemplateinfoline.o \
       readrestemplateureybradleyline.o \
       readrestemplatezmatline.o \
       readsoltemplateangleline.o \
       readsoltemplateatomline.o \
       readsoltemplatebondline.o \
       readsoltemplatecoordline.o \
       readsoltemplatedihedralline.o \
       readsoltemplateinfoline.o \
       readsoltemplateureybradleyline.o \
       readsoltemplatevariableline.o \
       readsolutetemplate.o \
       readsolventtemplate.o \
       readsvntemplateatomline.o \
       readsvntemplatecoordline.o \
       readsvntemplateinfoline.o \
       readtmpl.o \
       readub.o \
       recalculategeometries.o \
       recalculateparameters.o \
       resetgbenergy.o \
       residuebornintrapda.o \
       removeexcesscapsolvent.o \
       residuedetailedintraenergy.o \
       residueenergy.o \
       residuegcsolutenergy.o \
       residueintraenergy.o \
       residueintragbenergy.o \
       residueintranbenergy.o \
       residuemove.o \
       residueproteinenergy.o \
       residueresidueenergy.o \
       residueresidueborninterpda.o \
       residueresidueintergbenergy.o \
       residuesinglegcsoluteenergy.o \
       residuesinglesolventenergy.o \
       residuesolutebornpda.o \
       residuesoluteenergy.o \
       residuesolutegbenergy.o \
       residuesolventbornpda.o \
       residuesolventenergy.o \
       residuesolventgbenergy.o \
       residuesurface.o \
       retienergy.o \
       retitest.o	\
       reversevolumemove.o \
       revertborndata.o \
       revertcoordinates.o \
       revertenergies.o \
       rotamermove.o \
       rotateacid.o \
       rotategcsolute.o \
       rotatesolute.o \
       rotatesolvent.o \
       rotatevector.o \
       samplemove.o \
       savepoint.o \
       scanforvariablegroups.o \
       sepresresenergy.o \
       setallmoved.o \
       setbackwardsstructure.o \
       setdefaultoptions.o \
       setdirectangle.o \
       setdirectbond.o \
       setdirectdihedral.o \
       setdirectureybradley.o \
       setequilibriumgeometry.o \
       setforwardsstructure.o \
       setgcsolutemoved.o \
       setinitialtitrationstate.o \
       setkey.o \
       setlength.o \
       setproteingeometry.o \
       setresiduemoved.o \
       setsolutegeometry.o \
       setsolutemoved.o \
       setsolventmoved.o \
       setstream.o \
       setup.o \
       setupreplicas.o \
       setzmatrixgeometry.o \
       sgrnd.o \
       shigh.o \
       shutdown.o \
       sigma.o \
       simulate.o \
       simulation.o \
       singlepointenergy.o \
       savesolutegb.o \
       savesolventgb.o \
       saveresiduegb.o \
       slow.o \
       solutebornintrapda.o \
       solutedetailedintraenergy.o \
       soluteenergy.o \
       solutegcsenergy.o \
       soluteintraenergy.o \
       soluteintragbenergy.o \
       solutemove.o \
       solutesinglegcsoluteenergy.o \
       solutesinglesolventenergy.o \
       solutesoluteborninterpda.o \
       solutesolventbornpda.o \
       solutesoluteenergy.o \
       solutesoluteintergbenergy.o \
       solutesolventenergy.o \
       solutesolventgbenergy.o \
       solutesurface.o \
       solventbornintrapda.o \
       solventcapenergy.o \
       solventenergy.o \
       solventintragbenergy.o \
       solventsolventintergbenergy.o \
       solventmove.o \
       solventsolventborninterpda.o \
       solventsolventbulkenergy.o \
       solventsolventenergy.o \
       solventsolventsingleenergy.o \
       solventsurface.o \
       sortangle.o \
       sortbond.o \
       sortdihedral.o \
       sortmoltemplateangle.o \
       sortmoltemplatebond.o \
       sortmoltemplatedihedral.o \
       sortmoltemplateureybradley.o \
       sorttemplateangle.o \
       sorttemplatebond.o \
       sorttemplatedihedral.o \
       sorttemplateureybradley.o \
       sortureybradley.o \
       sortzangle.o \
       sortzbond.o \
       sortzdihedral.o \
       sortzureybradley.o \
       specifiedangle.o \
       specifiedbond.o \
       split.o \
       split2.o \
       splitchar.o \
       splitDDDSAsimulate.o \
       splitGBSAsimulate.o \
       stddev.o \
       storeborndata.o \
       storecoordinates.o \
       storeenergies.o \
       sumcharge.o \
       sumresiduecharge.o \
       sumsolutecharge.o \
       surfradius.o \
       surftension.o \
       switchbornradii.o \
       synchfepsolutes.o \
       syncproteinzmatrix.o \
       syncsolutedummies.o \
       takeacid.o \
       takebase.o \
       testenergies.o \
       teststring.o \
       thetamove.o \
       titrationmove.o \
       totalenergy.o \
       totalsurface.o \
       totalborn.o \
       translateacid.o \
       translategcsolute.o \
       translatesolute.o \
       translatesolvent.o \
       updateangle.o \
       updatebond.o \
       updatebornradii.o \
       updatedihedral.o \
       updateextraenergies.o \
       updategcsoluteenergies.o \
       updatemultiplemovesenergies.o \
       updateresidueenergies.o \
       updateresidueinternals.o \
       updatesoluteenergies.o \
       updatesoluteinternals.o \
       updatesolventenergies.o \
       updatesplitmovelist.o \
       updatesplitmovesolutelist.o \
       updatesplitmovesolventlist.o \
       updatesplitmoveresiduelist.o \
       updatesurfaces.o \
       updatetotalenergies.o \
       updateureybradley.o \
       ureybradleyenergy.o \
       veclength.o \
       veclength2.o \
       vecnormalise.o \
       volumemove.o \
       writenumberstring.o \
       writedoublestring.o	\
       writerestart.o \
       writezmatrix.o \
       zeroallpdacont.o \
       zeropdacont.o 

RETIOBJS = checkoutput.o	\
           fileexists.o	\
           fullsplit.o	\
           grnd.o	\
           getdateandtime.o	\
           getdouble.o	\
           getint.o	\
           getoption.o	\
           getretioptions.o	\
           initialiseprintsystem.o	\
           join.o	\
           killslaves.o	\
           linelength.o	\
           lowercase.o	\
           parsefile.o	\
           parseretiargs.o	\
           prettyline.o	\
           printline.o	\
           ranint.o	\
           readretienergy.o	\
           readspenergy.o	\
           retest.o	\
           reti_mpi.o	\
           retisimulation_mpi.o	\
           retiswap.o	\
           retitest.o	\
           sgrnd.o	\
           setkey.o	\
           spawnprocesses.o	\
           split.o	\
           tempswap.o	\
           writenumberstring.o

#
#*********************************************************************
#  You shouldn't need to modify anything below. These lines compile
#  ProtoMS based on the compiler names and flags, the object list
#*********************************************************************
#

# Dummy suffix for ftnchek targets:
.SUFFIXES: .chek
.PHONY: chekall

#escape sequence for a newline
nl = \033[\n

#
#
FC = $(fortrancompiler)
FC_MPI = $(fortrancompiler)

ifneq ($(HAVE_OPENMP),1)
  FFLAGS = $(fflags)
else
  FFLAGS = ${fflags} ${openmp_flags}
endif

# Flags for ftnchek:
CHEK1 = -pretty=no-long-line
CHEKFLAGS = $(CHEK1)

#colours...
# F/B light fg/bg; f/b dim fg/bg; Brgybpcw cols; 
# n return to normal; D dim; L light
n = \033[0m
D = \033[0m
L = \033[1m
FB = \033[1;30m
Fr = \033[1;31m
Fg = \033[1;32m
Fy = \033[1;33m
Fb = \033[1;34m
Fp = \033[1;35m
Fc = \033[1;36m
Fw = \033[1;37m
BB = \033[1;40m
Br = \033[1;41m
Bg = \033[1;42m
By = \033[1;43m
Bb = \033[1;44m
Bp = \033[1;45m
Bc = \033[1;46m
Bw = \033[1;47m
fB = \033[0;30m
fr = \033[0;31m
fg = \033[0;32m
fy = \033[0;33m
fb = \033[0;34m
fp = \033[0;35m
fc = \033[0;36m
fw = \033[0;37m
bB = \033[0;40m
br = \033[0;41m
bg = \033[0;42m
by = \033[0;43m
bb = \033[0;44m
bp = \033[0;45m 
bc = \033[0;46m 
bw = \033[0;47m

# Run ftnchek with all the .prj files except the one corresponding
# to the target's root:
%.chek : %.F ; \
  ftnchek $(filter-out $*.prj,$(PRJS)) $(CHEKFLAGS) \
    -noextern -library $<

# Derive a project file from a source file:
%.prj : %.F ; \
  ftnchek $(CHEKFLAGS) -noextern -project -library $<

# The list of objects is assumed to be in variable OBJS.
# Sources corresponding to the objects:
SRCS = $(OBJS:%.o=%.F)
# ftnchek project files:
PRJS = $(OBJS:%.o=%.prj)

ifneq ($(HAVE_MPI),0)
	fortrancompiler = $(mpi_fortrancompiler)
endif

%.o : %.F $(commons)
	@echo "$(Fp)Compiling $(Fb)$<...$(Fr)"
	@$(FC_MPI) -c $(FFLAGS) $< -o $@
	@echo "$(n)"

# Build the program
../protoms3: $(OBJS)
ifneq ($(NOPYTHON2),1)
	@echo "$(Fg)Generating compilation info...$(Fp)"
	@./getid.py idme.inc
endif
	@$(FC) -c $(FFLAGS) printheader.F -o printheader.o
	@echo "$(Fg)Linking the executable...$(Fr)"
	@$(FC) $(LFLAGS) $(FFLAGS) -o ../protoms3 $(OBJS) printheader.o

#ifneq ($(HAVE_MPI),0)
#	@echo "$(Fg)Making the RETI handler..."
#	@make ../reti
#endif

	@echo "$(Fg)Build complete. Go to the tutorials directory for"
	@echo "guidance on how to use this program.$(n)"


# Build the reti program
../reti : $(RETIOBJS)
	@echo "$(Fg)Linking the executable...$(Fr)"
	@$(FC_MPI) $(LFLAGS) $(FFLAGS) -o ../reti $(RETIOBJS)
	@echo "$(Fg)Build complete. Go to the tutorials directory for"
	@echo "guidance on how to use this program.$(n)"

chekall: $(PRJS) ; \
  ftnchek $(CHEKFLAGS) $(PRJS)

prjs: $(PRJS)

rebuild:
	@make distclean
	@make ../protoms3
	@make ../reti

#cleans everything!
clean:
	@echo -e "$(Fg)Making clean...$(Fp)"
	@-rm *.o
	@echo -e "$(Fc)Project now ready to be rebuilt or packaged...$(n)"

distclean:
	@make clean
	@-rm ../protoms3 ../reti

#rebuild the ftnchek files
cleancheck:
	@-rm *.prj
	@make chekall

# Rebuild dependencies:
depend: 
	@echo -e "$(Fg)Generating dependencies...$(Fr)"
	@./sfmakedepend $(SRCS)
	@echo -e "$(Fc)Please rebuild protomc via 'make rebuild'.$(n)"

# Make the documentation
docs:
	@echo -e "$(Fg)Generating autodocs...$(Fr)"
	@./parsefortran.py ../../manual/trunk/
	@echo -e "$(Fg)Generating manual.pdf...$(Fr)"
	@cd ../../manual/trunk/ ; make
	@echo -e "$(Fc)Manual now available as ../manual.pdf$(n)"
        
# DO NOT DELETE THIS LINE - used by make depend
addacid.o: proteins.inc dimensions.inc printer.inc
addacid.o:
addangle.o: dimensions.inc templates.inc
addangle.o:
addbase.o: proteins.inc dimensions.inc printer.inc
addbase.o:
addbond.o: dimensions.inc templates.inc
addbond.o:
addbondharmonicrestraint.o: dimensions.inc extraenergies.inc gb.inc surface.inc
addbondharmonicrestraint.o: constants.inc
addbondharmonicrestraint.o:
addcartesianflatbottomrestraint.o: dimensions.inc extraenergies.inc gb.inc
addcartesianflatbottomrestraint.o: surface.inc constants.inc
addcartesianflatbottomrestraint.o:
addcartesianharmonicrestraint.o: dimensions.inc extraenergies.inc gb.inc
addcartesianharmonicrestraint.o: surface.inc constants.inc
addcartesianharmonicrestraint.o:
addchainangle.o: dimensions.inc templates.inc
addchainangle.o:
addchainbond.o: dimensions.inc templates.inc
addchainbond.o:
addchaindihedral.o: enums.inc proteins.inc dimensions.inc printer.inc
addchaindihedral.o: templates.inc
addchaindihedral.o:
addchainureybradley.o: dimensions.inc templates.inc
addchainureybradley.o:
addclj.o: dimensions.inc constants.inc parameters.inc printer.inc
addclj.o:
adddihedral.o: enums.inc proteins.inc dimensions.inc printer.inc templates.inc
adddihedral.o:
adddihedralharmonicrestraint.o: dimensions.inc extraenergies.inc gb.inc
adddihedralharmonicrestraint.o: surface.inc constants.inc printer.inc
adddihedralharmonicrestraint.o:
addid.o: enums.inc proteins.inc dimensions.inc solutes.inc solvents.inc
addid.o: printer.inc splitter.inc labels.inc
addid.o:
addmolangle.o: enums.inc proteins.inc dimensions.inc printer.inc templates.inc
addmolangle.o:
addmolatom.o: enums.inc dimensions.inc solutes.inc printer.inc templates.inc
addmolatom.o:
addmolbond.o: enums.inc proteins.inc dimensions.inc printer.inc templates.inc
addmolbond.o:
addmoldihedral.o: enums.inc proteins.inc dimensions.inc printer.inc
addmoldihedral.o: templates.inc
addmoldihedral.o:
addmolureybradley.o: enums.inc dimensions.inc printer.inc templates.inc
addmolureybradley.o:
addproteinangle.o: proteins.inc enums.inc dimensions.inc printer.inc
addproteinangle.o:
addproteinbond.o: proteins.inc enums.inc dimensions.inc printer.inc
addproteinbond.o:
addproteindihedral.o: proteins.inc enums.inc dimensions.inc printer.inc
addproteindihedral.o:
addproteinureybradley.o: proteins.inc enums.inc dimensions.inc printer.inc
addproteinureybradley.o:
addsoluteangle.o: enums.inc dimensions.inc solutes.inc printer.inc
addsoluteangle.o:
addsolutebond.o: enums.inc dimensions.inc solutes.inc printer.inc
addsolutebond.o:
addsolutedihedral.o: enums.inc dimensions.inc solutes.inc printer.inc
addsolutedihedral.o:
addsolutedummies.o: enums.inc dimensions.inc solutes.inc
addsolutedummies.o:
addsolutenb.o: enums.inc dimensions.inc solutes.inc printer.inc
addsolutenb.o:
addsoluteureybradley.o: enums.inc dimensions.inc solutes.inc printer.inc
addsoluteureybradley.o:
addsolventaxes.o: simparams.inc enums.inc dimensions.inc solvents.inc
addsolventaxes.o:
addureybradley.o: dimensions.inc templates.inc
addureybradley.o:
amber.o: simparams.inc dimensions.inc parameters.inc
amber.o:
angleenergy.o: enums.inc proteins.inc dimensions.inc solutes.inc parameters.inc
angleenergy.o: constants.inc printer.inc flags.inc
angleenergy.o:
assignborn.o: dimensions.inc gb.inc constants.inc parameters.inc printer.inc
assignborn.o:
assignresidueoverrides.o: proteins.inc dimensions.inc printer.inc templates.inc
assignresidueoverrides.o:
assignresiduetemplate.o: enums.inc proteins.inc dimensions.inc constants.inc
assignresiduetemplate.o: printer.inc templates.inc
assignresiduetemplate.o:
assignrotamerlibrary.o: proteins.inc dimensions.inc printer.inc rotamer.inc
assignrotamerlibrary.o:
assignsolutetemplates.o: enums.inc simparams.inc dimensions.inc solutes.inc
assignsolutetemplates.o: constants.inc printer.inc templates.inc
assignsolutetemplates.o:
assignsolutevariables.o: simparams.inc dimensions.inc solutes.inc printer.inc
assignsolutevariables.o: templates.inc
assignsolutevariables.o:
assignsolventtemplate.o: enums.inc dimensions.inc solvents.inc printer.inc
assignsolventtemplate.o: templates.inc
assignsolventtemplate.o:
assignsolventtemplates.o: dimensions.inc solvents.inc
assignsolventtemplates.o:
assignsurface.o: dimensions.inc constants.inc parameters.inc surface.inc
assignsurface.o: printer.inc
assignsurface.o:
assignsystem.o: enums.inc proteins.inc dimensions.inc solutes.inc constants.inc
assignsystem.o: solvents.inc printer.inc simparams.inc flags.inc pressure.inc
assignsystem.o:
atomname.o: proteins.inc dimensions.inc
atomname.o:
bndharmonicrestraintenergy.o: simparams.inc flags.inc enums.inc dimensions.inc
bndharmonicrestraintenergy.o: extraenergies.inc gb.inc surface.inc constants.inc
bndharmonicrestraintenergy.o: printer.inc movelist.inc
bndharmonicrestraintenergy.o:
bondedresresenergy.o: nbenergy2.inc proteins.inc enums.inc dimensions.inc
bondedresresenergy.o: parameters.inc constants.inc nbpropro.inc printer.inc
bondedresresenergy.o: nbenergyvar2.inc simparams.inc flags.inc quickpot.inc
bondedresresenergy.o: nbinitial.inc
bondedresresenergy.o:
bondenergy.o: enums.inc proteins.inc dimensions.inc solutes.inc parameters.inc
bondenergy.o: constants.inc printer.inc flags.inc
bondenergy.o:
buildprotein.o: proteins.inc enums.inc dimensions.inc printer.inc
buildprotein.o:
buildresidue.o: proteins.inc dimensions.inc
buildresidue.o:
buildsolute.o: enums.inc dimensions.inc solutes.inc
buildsolute.o:
buildsolventbox.o: simparams.inc proteins.inc dimensions.inc solutes.inc
buildsolventbox.o: constants.inc solvents.inc printer.inc
buildsolventbox.o:
calcangle.o: enums.inc simparams.inc proteins.inc dimensions.inc solutes.inc
calcangle.o: constants.inc printer.inc
calcangle.o:
calcaverageenergy.o: proteins.inc dimensions.inc solutes.inc solvents.inc
calcaverageenergy.o: constants.inc printer.inc averages.inc flags.inc
calcaverageenergy.o: energies.inc metropolis.inc extraenergies.inc gb.inc
calcaverageenergy.o: surface.inc pressure.inc
calcaverageenergy.o:
calcaverageextraenergies.o: flags.inc dimensions.inc extraenergies.inc gb.inc
calcaverageextraenergies.o: surface.inc constants.inc printer.inc
calcaverageextraenergies.o:
calcaveragegbenergy.o: dimensions.inc gb.inc constants.inc
calcaveragegbenergy.o:
calcaveragesaenergy.o: dimensions.inc constants.inc surface.inc
calcaveragesaenergy.o:
calcbond.o: enums.inc simparams.inc proteins.inc dimensions.inc solutes.inc
calcbond.o: constants.inc printer.inc
calcbond.o:
calcdihedral.o: proteins.inc enums.inc dimensions.inc solutes.inc constants.inc
calcdihedral.o: printer.inc
calcdihedral.o:
calcgeneraldihedral.o: enums.inc dimensions.inc constants.inc printer.inc
calcgeneraldihedral.o:
calculatecljpairs.o: flags.inc dimensions.inc constants.inc parameters.inc
calculatecljpairs.o: printer.inc
calculatecljpairs.o:
calculatesolventbox.o: simparams.inc flags.inc dimensions.inc solvents.inc
calculatesolventbox.o:
calculatetotalsurface.o: enums.inc proteins.inc dimensions.inc solutes.inc
calculatetotalsurface.o: solvents.inc constants.inc printer.inc surface.inc
calculatetotalsurface.o: simparams.inc
calculatetotalsurface.o:
calcureybradley.o: enums.inc simparams.inc proteins.inc dimensions.inc
calcureybradley.o: solutes.inc constants.inc printer.inc
calcureybradley.o:
changelambda.o: simparams.inc flags.inc enums.inc energies.inc dimensions.inc
changelambda.o: printer.inc
changelambda.o:
charge.o: dimensions.inc constants.inc parameters.inc
charge.o:
checkbonded.o: proteins.inc dimensions.inc constants.inc
checkbonded.o:
checkgottotal.o: enums.inc flags.inc dimensions.inc gb.inc printer.inc
checkgottotal.o:
cleansoluteparameters.o: dimensions.inc solutes.inc constants.inc parameters.inc
cleansoluteparameters.o:
clearid.o: dimensions.inc splitter.inc labels.inc printer.inc
clearid.o:
clearmovelist.o: dimensions.inc movelist.inc
clearmovelist.o:
clearresiduecoordstack.o: dimensions.inc stacks.inc
clearresiduecoordstack.o:
clearresidueenergystack.o: enums.inc dimensions.inc printer.inc stacks.inc
clearresidueenergystack.o:
clearsaveresiduegb.o: dimensions.inc gb.inc
clearsaveresiduegb.o:
clearsavesolutegb.o: dimensions.inc gb.inc
clearsavesolutegb.o:
clearsavesolventgb.o: dimensions.inc gb.inc
clearsavesolventgb.o:
clearsolutecoordstack.o: dimensions.inc stacks.inc
clearsolutecoordstack.o:
clearsolventcoordstack.o: dimensions.inc stacks.inc
clearsolventcoordstack.o:
clearsplitmovelist.o: dimensions.inc quickpot.inc
clearsplitmovelist.o:
clearsurfaces.o: proteins.inc dimensions.inc solutes.inc solvents.inc
clearsurfaces.o: constants.inc surface.inc printer.inc simparams.inc
clearsurfaces.o:
cljid.o: simparams.inc dimensions.inc constants.inc parameters.inc
cljid.o:
closefile.o: dimensions.inc printer.inc
closefile.o:
closeprogram.o: printer.inc
closeprogram.o:
convertres.o: enums.inc proteins.inc dimensions.inc constants.inc printer.inc
convertres.o: templates.inc
convertres.o:
copysolvent.o: dimensions.inc solvents.inc printer.inc
copysolvent.o:
createconnectivitymap.o: proteins.inc flags.inc enums.inc dimensions.inc
createconnectivitymap.o: printer.inc
createconnectivitymap.o:
createresresconnectivitymap.o: proteins.inc flags.inc dimensions.inc printer.inc
createresresconnectivitymap.o:
deprotonate.o: proteins.inc dimensions.inc
deprotonate.o:
detailedresidueenergy.o: simparams.inc proteins.inc enums.inc dimensions.inc
detailedresidueenergy.o: solutes.inc solvents.inc constants.inc printer.inc
detailedresidueenergy.o:
detailedsoluteenergy.o: enums.inc proteins.inc dimensions.inc solutes.inc
detailedsoluteenergy.o: solvents.inc constants.inc printer.inc
detailedsoluteenergy.o:
dihedralenergy.o: enums.inc proteins.inc dimensions.inc solutes.inc
dihedralenergy.o: parameters.inc constants.inc printer.inc simparams.inc
dihedralenergy.o: flags.inc
dihedralenergy.o:
dihedralharmonicrestraintenergy.o: enums.inc proteins.inc dimensions.inc
dihedralharmonicrestraintenergy.o: solutes.inc constants.inc solvents.inc
dihedralharmonicrestraintenergy.o: printer.inc flags.inc extraenergies.inc
dihedralharmonicrestraintenergy.o: gb.inc surface.inc movelist.inc
dihedralharmonicrestraintenergy.o:
dihenergy.o: dimensions.inc constants.inc parameters.inc
dihenergy.o:
dualtopologysasa.o: proteins.inc enums.inc dimensions.inc solutes.inc
dualtopologysasa.o: solvents.inc constants.inc surface.inc printer.inc flags.inc
dualtopologysasa.o: simparams.inc energies.inc extraenergies.inc gb.inc
dualtopologysasa.o: surface.inc
dualtopologysasa.o:
duplicatesolute.o: dimensions.inc solutes.inc printer.inc
duplicatesolute.o:
epslon.o: dimensions.inc constants.inc parameters.inc
epslon.o:
equilibrate.o: enums.inc dimensions.inc printer.inc
equilibrate.o:
extraenergy.o: enums.inc flags.inc dimensions.inc extraenergies.inc gb.inc
extraenergy.o: surface.inc solutes.inc constants.inc printer.inc
extraenergy.o:
fakesimulate.o: enums.inc
fakesimulate.o:
findSA.o: dimensions.inc constants.inc surface.inc
findSA.o:
findangleindex.o: dimensions.inc parameters.inc
findangleindex.o:
findangleparam.o: enums.inc proteins.inc dimensions.inc parameters.inc
findangleparam.o: printer.inc
findangleparam.o:
findatom.o: proteins.inc dimensions.inc
findatom.o:
findbondindex.o: dimensions.inc parameters.inc
findbondindex.o:
findbondparam.o: enums.inc proteins.inc dimensions.inc parameters.inc
findbondparam.o: printer.inc
findbondparam.o:
findchainangle.o: dimensions.inc templates.inc
findchainangle.o:
findchainatom.o: enums.inc dimensions.inc printer.inc templates.inc
findchainatom.o:
findchainbond.o: dimensions.inc templates.inc
findchainbond.o:
findchaindihedral.o: dimensions.inc templates.inc
findchaindihedral.o:
findchaintemplate.o: dimensions.inc templates.inc
findchaintemplate.o:
findchainureybradley.o: dimensions.inc templates.inc
findchainureybradley.o:
finddihedralindex.o: dimensions.inc parameters.inc printer.inc
finddihedralindex.o:
finddihedralparam.o: enums.inc proteins.inc dimensions.inc parameters.inc
finddihedralparam.o: printer.inc
finddihedralparam.o:
findgb.o: dimensions.inc gb.inc constants.inc
findgb.o:
findmoleculetemplate.o: dimensions.inc printer.inc templates.inc
findmoleculetemplate.o:
findmoltemplateangle.o: dimensions.inc templates.inc
findmoltemplateangle.o:
findmoltemplateatom.o: dimensions.inc printer.inc templates.inc
findmoltemplateatom.o:
findmoltemplatebond.o: dimensions.inc templates.inc
findmoltemplatebond.o:
findmoltemplatedihedral.o: dimensions.inc templates.inc
findmoltemplatedihedral.o:
findmoltemplateureybradley.o: dimensions.inc templates.inc
findmoltemplateureybradley.o:
findoverride.o: dimensions.inc templates.inc
findoverride.o:
findresidue.o: proteins.inc dimensions.inc
findresidue.o:
findresiduetemplate.o: dimensions.inc templates.inc
findresiduetemplate.o:
findsolresidue.o: enums.inc dimensions.inc solutes.inc printer.inc
findsolresidue.o:
findsoluteatom.o: enums.inc dimensions.inc solutes.inc printer.inc
findsoluteatom.o:
findsolventatom.o: dimensions.inc solvents.inc
findsolventatom.o:
findsolventtemplate.o: dimensions.inc templates.inc
findsolventtemplate.o:
findsolventtemplateatom.o: dimensions.inc templates.inc
findsolventtemplateatom.o:
findtemplate.o: dimensions.inc templates.inc
findtemplate.o:
findtemplateangle.o: dimensions.inc templates.inc
findtemplateangle.o:
findtemplateatom.o: enums.inc dimensions.inc printer.inc templates.inc
findtemplateatom.o:
findtemplatebond.o: dimensions.inc templates.inc
findtemplatebond.o:
findtemplatedihedral.o: dimensions.inc templates.inc
findtemplatedihedral.o:
findtemplateureybradley.o: dimensions.inc templates.inc
findtemplateureybradley.o:
findureybradleyindex.o: dimensions.inc parameters.inc
findureybradleyindex.o:
findureybradleyparam.o: enums.inc dimensions.inc parameters.inc printer.inc
findureybradleyparam.o:
flagmovedbits.o: dimensions.inc quickpot.inc solutes.inc movelist.inc
flagmovedbits.o:
flatbottomrestraintenergy.o: enums.inc dimensions.inc extraenergies.inc gb.inc
flatbottomrestraintenergy.o: surface.inc constants.inc printer.inc movelist.inc
flatbottomrestraintenergy.o:
fliphistidine.o: proteins.inc dimensions.inc printer.inc
fliphistidine.o:
gbenergy.o: enums.inc proteins.inc dimensions.inc solutes.inc gb.inc
gbenergy.o: solvents.inc constants.inc printer.inc flags.inc simparams.inc
gbenergy.o: quickpot.inc
gbenergy.o:
gbsavesystem.o: enums.inc proteins.inc dimensions.inc solutes.inc gb.inc
gbsavesystem.o: solvents.inc constants.inc printer.inc flags.inc simparams.inc
gbsavesystem.o:
genproteinatom.o: proteins.inc dimensions.inc
genproteinatom.o:
gensoluteatom.o: dimensions.inc solutes.inc
gensoluteatom.o:
getangle.o: dimensions.inc parameters.inc printer.inc
getangle.o:
getatomname.o: proteins.inc dimensions.inc
getatomname.o:
getaxisvectors.o: dimensions.inc constants.inc printer.inc
getaxisvectors.o:
getbond.o: dimensions.inc parameters.inc printer.inc
getbond.o:
getchaintemplate.o: enums.inc dimensions.inc templates.inc
getchaintemplate.o:
getclj.o: simparams.inc dimensions.inc gb.inc constants.inc parameters.inc
getclj.o: surface.inc printer.inc
getclj.o:
getcoordinates.o: proteins.inc enums.inc dimensions.inc solutes.inc solvents.inc
getcoordinates.o: printer.inc
getcoordinates.o:
getdihedral.o: dimensions.inc parameters.inc printer.inc
getdihedral.o:
getdihwells.o: dimensions.inc splitter.inc printer.inc
getdihwells.o:
getdofinfo.o: dimensions.inc splitter.inc
getdofinfo.o:
geteffectivebornradii.o: proteins.inc dimensions.inc solutes.inc gb.inc
geteffectivebornradii.o: solvents.inc parameters.inc constants.inc printer.inc
geteffectivebornradii.o: flags.inc
geteffectivebornradii.o:
getgbiborn.o: dimensions.inc gb.inc constants.inc parameters.inc printer.inc
getgbiborn.o:
getgbscale.o: dimensions.inc gb.inc constants.inc parameters.inc printer.inc
getgbscale.o:
getmajorminoraxis.o: proteins.inc enums.inc dimensions.inc solutes.inc
getmajorminoraxis.o: solvents.inc constants.inc printer.inc
getmajorminoraxis.o:
getmoveinfo.o: enums.inc proteins.inc dimensions.inc solutes.inc solvents.inc
getmoveinfo.o: constants.inc printer.inc moves.inc simparams.inc splitter.inc
getmoveinfo.o: pressure.inc
getmoveinfo.o:
getmovetype.o: enums.inc dimensions.inc printer.inc moves.inc
getmovetype.o:
getnonbonded.o: enums.inc proteins.inc dimensions.inc printer.inc
getnonbonded.o:
getoption.o: dimensions.inc splitter.inc keyset.inc
getoption.o:
getoptions.o: enums.inc proteins.inc dimensions.inc solutes.inc gb.inc
getoptions.o: parameters.inc solvents.inc constants.inc surface.inc printer.inc
getoptions.o: flags.inc simparams.inc metropolis.inc splitter.inc pressure.inc
getoptions.o:
getperpendicular.o: dimensions.inc constants.inc
getperpendicular.o:
getresiduecentreradii.o: proteins.inc enums.inc dimensions.inc printer.inc
getresiduecentreradii.o:
getresiduecentreradius.o: proteins.inc dimensions.inc constants.inc
getresiduecentreradius.o:
getresintranonbonded.o: proteins.inc flags.inc dimensions.inc printer.inc
getresintranonbonded.o:
getresresbondindex.o: proteins.inc dimensions.inc
getresresbondindex.o:
getresresinternonbonded.o: proteins.inc flags.inc dimensions.inc printer.inc
getresresinternonbonded.o:
getsmallestdelta.o: dimensions.inc constants.inc printer.inc
getsmallestdelta.o:
getsolatomname.o: dimensions.inc solutes.inc
getsolatomname.o:
getsolutecentreradius.o: dimensions.inc solutes.inc constants.inc
getsolutecentreradius.o:
getsolutenonbonded.o: enums.inc dimensions.inc solutes.inc constants.inc
getsolutenonbonded.o: printer.inc
getsolutenonbonded.o:
getsolutezmatrix.o: enums.inc dimensions.inc solutes.inc
getsolutezmatrix.o:
getureybradley.o: dimensions.inc parameters.inc
getureybradley.o:
getzmatdihedral.o: proteins.inc enums.inc flags.inc dimensions.inc solutes.inc
getzmatdihedral.o: solvents.inc printer.inc
getzmatdihedral.o:
harmonicrestraintenergy.o: enums.inc dimensions.inc extraenergies.inc gb.inc
harmonicrestraintenergy.o: surface.inc constants.inc printer.inc movelist.inc
harmonicrestraintenergy.o:
hasangle.o: proteins.inc dimensions.inc
hasangle.o:
hasatommoved.o: enums.inc proteins.inc dimensions.inc solutes.inc solvents.inc
hasatommoved.o: movelist.inc
hasatommoved.o:
hasbond.o: proteins.inc dimensions.inc
hasbond.o:
hasdihedral.o: proteins.inc dimensions.inc
hasdihedral.o:
hasresiduemoved.o: movelist.inc
hasresiduemoved.o:
hassoluteangle.o: dimensions.inc solutes.inc
hassoluteangle.o:
hassolutebond.o: dimensions.inc solutes.inc
hassolutebond.o:
hassolutedihedral.o: dimensions.inc solutes.inc
hassolutedihedral.o:
hassolutemoved.o: movelist.inc
hassolutemoved.o:
hassoluteureybradley.o: dimensions.inc solutes.inc
hassoluteureybradley.o:
hassolventmoved.o: movelist.inc
hassolventmoved.o:
hasureybradley.o: proteins.inc dimensions.inc
hasureybradley.o:
initialise.o: version.inc dimensions.inc printer.inc
initialise.o:
initialisearrays.o: proteins.inc dimensions.inc solutes.inc solvents.inc
initialisearrays.o: constants.inc moves.inc quickpot.inc labels.inc
initialisearrays.o: shrakeinit.inc enums.inc parameters.inc rotamer.inc
initialisearrays.o: printer.inc averages.inc flags.inc simparams.inc
initialisearrays.o: metropolis.inc energies.inc extraenergies.inc gb.inc
initialisearrays.o: surface.inc pressure.inc templates.inc stacks.inc
initialisearrays.o:
initialiseextraenergies.o: dimensions.inc extraenergies.inc gb.inc surface.inc
initialiseextraenergies.o: constants.inc
initialiseextraenergies.o:
initialiseprintsystem.o: dimensions.inc printer.inc
initialiseprintsystem.o:
initialiseprotein.o: proteins.inc dimensions.inc constants.inc
initialiseprotein.o:
initialisesolute.o: dimensions.inc solutes.inc constants.inc moves.inc
initialisesolute.o:
inresiduesplitmoved.o: dimensions.inc quickpot.inc printer.inc
inresiduesplitmoved.o:
insolutesplitmoved.o: dimensions.inc quickpot.inc printer.inc
insolutesplitmoved.o:
insolventsplitmoved.o: dimensions.inc quickpot.inc printer.inc
insolventsplitmoved.o:
interpretposition.o: enums.inc
interpretposition.o:
interresidueresidueenergy.o: nbenergy2.inc proteins.inc enums.inc dimensions.inc
interresidueresidueenergy.o: parameters.inc constants.inc nbpropro.inc
interresidueresidueenergy.o: printer.inc nbenergyvar2.inc simparams.inc
interresidueresidueenergy.o: flags.inc quickpot.inc nbinitial.inc
interresidueresidueenergy.o:
intraresidueproteinenergy.o: proteins.inc dimensions.inc constants.inc
intraresidueproteinenergy.o:
isbackbone.o: proteins.inc dimensions.inc
isbackbone.o:
join.o: dimensions.inc splitter.inc
join.o:
lambdamove.o: dimensions.inc printer.inc
lambdamove.o:
loadprotein.o: enums.inc proteins.inc dimensions.inc splitter.inc constants.inc
loadprotein.o: printer.inc
loadprotein.o:
loadsolute.o: enums.inc dimensions.inc splitter.inc solutes.inc constants.inc
loadsolute.o: printer.inc
loadsolute.o:
loadsolvent.o: enums.inc dimensions.inc solvents.inc constants.inc printer.inc
loadsolvent.o: simparams.inc flags.inc splitter.inc
loadsolvent.o:
main.o: protoms.inc proteins.inc dimensions.inc solutes.inc solvents.inc
main.o: constants.inc keyset.inc moves.inc movelist.inc enums.inc parameters.inc
main.o: printer.inc averages.inc simparams.inc flags.inc metropolis.inc
main.o: energies.inc splitter.inc extraenergies.inc gb.inc surface.inc
main.o: pressure.inc stacks.inc templates.inc
main.o:
makevariableangle.o: proteins.inc dimensions.inc printer.inc
makevariableangle.o:
makevariablebond.o: proteins.inc dimensions.inc printer.inc
makevariablebond.o:
makevariabledihedral.o: proteins.inc dimensions.inc printer.inc
makevariabledihedral.o:
makevariablesolangle.o: enums.inc proteins.inc dimensions.inc solutes.inc
makevariablesolangle.o: printer.inc
makevariablesolangle.o:
makevariablesolbond.o: enums.inc proteins.inc dimensions.inc solutes.inc
makevariablesolbond.o: printer.inc
makevariablesolbond.o:
makevariablesoldihedral.o: enums.inc proteins.inc dimensions.inc solutes.inc
makevariablesoldihedral.o: printer.inc
makevariablesoldihedral.o:
markgbupdate.o: enums.inc proteins.inc dimensions.inc solutes.inc gb.inc
markgbupdate.o: solvents.inc constants.inc printer.inc flags.inc simparams.inc
markgbupdate.o:
marksaupdate.o: proteins.inc dimensions.inc solutes.inc constants.inc
marksaupdate.o: solvents.inc surface.inc printer.inc flags.inc energies.inc
marksaupdate.o:
metropolis.o: enums.inc proteins.inc dimensions.inc solutes.inc solvents.inc
metropolis.o: constants.inc printer.inc metropolis.inc pressure.inc
metropolis.o:
minsolventproteindist2.o: proteins.inc flags.inc simparams.inc enums.inc
minsolventproteindist2.o: dimensions.inc solvents.inc
minsolventproteindist2.o:
minsolventsolutedist2.o: flags.inc simparams.inc enums.inc dimensions.inc
minsolventsolutedist2.o: solutes.inc solvents.inc
minsolventsolutedist2.o:
moveacid.o: proteins.inc enums.inc dimensions.inc constants.inc printer.inc
moveacid.o:
moverotamer.o: proteins.inc dimensions.inc constants.inc printer.inc rotamer.inc
moverotamer.o:
movesidechain.o: proteins.inc dimensions.inc
movesidechain.o:
movesolute.o: dimensions.inc solutes.inc constants.inc printer.inc
movesolute.o:
movesolvent.o: flags.inc simparams.inc enums.inc dimensions.inc constants.inc
movesolvent.o: solvents.inc
movesolvent.o:
movevolume.o: enums.inc proteins.inc dimensions.inc solutes.inc solvents.inc
movevolume.o: constants.inc flags.inc simparams.inc pressure.inc
movevolume.o:
multiplemovesenergy.o: enums.inc proteins.inc dimensions.inc solutes.inc
multiplemovesenergy.o: constants.inc solvents.inc printer.inc flags.inc
multiplemovesenergy.o: simparams.inc energies.inc movelist.inc
multiplemovesenergy.o:
namestring.o: proteins.inc dimensions.inc
namestring.o:
notbonded.o: proteins.inc enums.inc dimensions.inc solutes.inc
notbonded.o:
openfile.o: dimensions.inc printer.inc
openfile.o:
openstreamfiles.o: dimensions.inc splitter.inc printer.inc
openstreamfiles.o:
parseargs.o: version.inc dimensions.inc splitter.inc printer.inc
parseargs.o:
parsefile.o: dimensions.inc splitter.inc printer.inc keyset.inc
parsefile.o:
perturbforcefield.o: dimensions.inc constants.inc
perturbforcefield.o:
popallresidueenergy.o: enums.inc proteins.inc energies.inc dimensions.inc
popallresidueenergy.o: solutes.inc printer.inc stacks.inc
popallresidueenergy.o:
popdeletesolventcoords.o: enums.inc dimensions.inc printer.inc stacks.inc
popdeletesolventcoords.o:
poppdacont.o: proteins.inc dimensions.inc solutes.inc gb.inc solvents.inc
poppdacont.o: parameters.inc constants.inc
poppdacont.o:
popresiduecoords.o: proteins.inc dimensions.inc printer.inc stacks.inc
popresiduecoords.o:
popsolutecoords.o: enums.inc dimensions.inc solutes.inc printer.inc stacks.inc
popsolutecoords.o:
popsolventcoords.o: enums.inc dimensions.inc solvents.inc printer.inc stacks.inc
popsolventcoords.o:
positionstring.o: enums.inc
positionstring.o:
prettyline.o: dimensions.inc printer.inc
prettyline.o:
printattemptaccept.o: enums.inc proteins.inc dimensions.inc solutes.inc
printattemptaccept.o: solvents.inc printer.inc moves.inc pressure.inc
printattemptaccept.o:
printaverageextraenergies.o: flags.inc dimensions.inc extraenergies.inc gb.inc
printaverageextraenergies.o: surface.inc constants.inc printer.inc
printaverageextraenergies.o:
printaverages.o: enums.inc proteins.inc dimensions.inc solutes.inc constants.inc
printaverages.o: solvents.inc parameters.inc printer.inc averages.inc flags.inc
printaverages.o: simparams.inc metropolis.inc pressure.inc
printaverages.o:
printaveragesaenergy.o: dimensions.inc constants.inc surface.inc printer.inc
printaveragesaenergy.o:
printaveragesgbenergy.o: dimensions.inc gb.inc constants.inc printer.inc
printaveragesgbenergy.o:
printboxtext.o: dimensions.inc splitter.inc
printboxtext.o:
printenergies.o: energies.inc dimensions.inc printer.inc
printenergies.o:
printline.o: dimensions.inc printer.inc
printline.o:
printparameters.o: simparams.inc proteins.inc dimensions.inc solutes.inc
printparameters.o: parameters.inc solvents.inc printer.inc
printparameters.o:
printpdb.o: enums.inc proteins.inc dimensions.inc solutes.inc solvents.inc
printpdb.o: constants.inc printer.inc simparams.inc flags.inc
printpdb.o:
printpdbatom.o: dimensions.inc parameters.inc printer.inc
printpdbatom.o:
printpdbprotein.o: proteins.inc enums.inc dimensions.inc printer.inc
printpdbprotein.o:
printpdbresidue.o: proteins.inc dimensions.inc
printpdbresidue.o:
printpdbsolute.o: dimensions.inc solutes.inc printer.inc
printpdbsolute.o:
printpdbsolvent.o: dimensions.inc solvents.inc printer.inc
printpdbsolvent.o:
processgeometryline.o: enums.inc proteins.inc dimensions.inc solutes.inc
processgeometryline.o: constants.inc solvents.inc printer.inc splitter.inc
processgeometryline.o:
processidline.o: dimensions.inc splitter.inc printer.inc
processidline.o:
processrestraintline.o: enums.inc proteins.inc dimensions.inc solutes.inc
processrestraintline.o: constants.inc solvents.inc printer.inc labels.inc
processrestraintline.o: splitter.inc
processrestraintline.o:
processsolventline.o: dimensions.inc splitter.inc constants.inc printer.inc
processsolventline.o:
proteinintraenergy.o: proteins.inc enums.inc dimensions.inc constants.inc
proteinintraenergy.o:
proteinproteinenergy.o: proteins.inc enums.inc dimensions.inc constants.inc
proteinproteinenergy.o:
proteinresresbondlist.o: proteins.inc dimensions.inc printer.inc
proteinresresbondlist.o:
proteinsoluteenergy.o: enums.inc flags.inc simparams.inc proteins.inc
proteinsoluteenergy.o: dimensions.inc solutes.inc constants.inc
proteinsoluteenergy.o:
proteinsolventenergy.o: proteins.inc enums.inc dimensions.inc constants.inc
proteinsolventenergy.o: solvents.inc
proteinsolventenergy.o:
protonate.o: proteins.inc dimensions.inc
protonate.o:
protonnum.o: simparams.inc dimensions.inc parameters.inc
protonnum.o:
pushbacknewresidueenergy.o: enums.inc proteins.inc energies.inc dimensions.inc
pushbacknewresidueenergy.o: solutes.inc printer.inc stacks.inc
pushbacknewresidueenergy.o:
pushbackoldresidueenergy.o: enums.inc proteins.inc energies.inc dimensions.inc
pushbackoldresidueenergy.o: solutes.inc printer.inc stacks.inc
pushbackoldresidueenergy.o:
pushbackpdacont.o: proteins.inc dimensions.inc solutes.inc gb.inc solvents.inc
pushbackpdacont.o: parameters.inc constants.inc
pushbackpdacont.o:
pushbackresiduecoords.o: proteins.inc dimensions.inc printer.inc stacks.inc
pushbackresiduecoords.o:
pushbacksolutecoords.o: enums.inc dimensions.inc solutes.inc printer.inc
pushbacksolutecoords.o: stacks.inc
pushbacksolutecoords.o:
pushbacksolventcoords.o: enums.inc dimensions.inc solvents.inc printer.inc
pushbacksolventcoords.o: stacks.inc
pushbacksolventcoords.o:
readang.o: enums.inc proteins.inc dimensions.inc splitter.inc parameters.inc
readang.o: constants.inc printer.inc
readang.o:
readbnd.o: enums.inc proteins.inc dimensions.inc splitter.inc parameters.inc
readbnd.o: constants.inc printer.inc
readbnd.o:
readborn.o: dimensions.inc splitter.inc gb.inc constants.inc printer.inc
readborn.o:
readchaintemplate.o: enums.inc dimensions.inc splitter.inc constants.inc
readchaintemplate.o: printer.inc templates.inc
readchaintemplate.o:
readchaintemplateangleline.o: enums.inc proteins.inc dimensions.inc splitter.inc
readchaintemplateangleline.o: printer.inc templates.inc
readchaintemplateangleline.o:
readchaintemplateatomline.o: enums.inc dimensions.inc splitter.inc printer.inc
readchaintemplateatomline.o: templates.inc
readchaintemplateatomline.o:
readchaintemplatebbatomline.o: dimensions.inc splitter.inc printer.inc
readchaintemplatebbatomline.o: templates.inc
readchaintemplatebbatomline.o:
readchaintemplatebondline.o: enums.inc proteins.inc dimensions.inc splitter.inc
readchaintemplatebondline.o: printer.inc templates.inc
readchaintemplatebondline.o:
readchaintemplatedihedralline.o: enums.inc proteins.inc dimensions.inc
readchaintemplatedihedralline.o: splitter.inc printer.inc templates.inc
readchaintemplatedihedralline.o:
readchaintemplateubline.o: dimensions.inc splitter.inc printer.inc templates.inc
readchaintemplateubline.o:
readchaintemplatezmatline.o: dimensions.inc splitter.inc printer.inc
readchaintemplatezmatline.o: templates.inc
readchaintemplatezmatline.o:
readclj.o: enums.inc dimensions.inc splitter.inc parameters.inc constants.inc
readclj.o: printer.inc
readclj.o:
readdih.o: enums.inc proteins.inc dimensions.inc splitter.inc parameters.inc
readdih.o: constants.inc printer.inc
readdih.o:
readdouble.o: dimensions.inc printer.inc
readdouble.o:
readinfo.o: flags.inc simparams.inc dimensions.inc splitter.inc constants.inc
readinfo.o: printer.inc
readinfo.o:
readinfoline.o: dimensions.inc splitter.inc
readinfoline.o:
readint.o: dimensions.inc printer.inc
readint.o:
readlogical.o: dimensions.inc printer.inc
readlogical.o:
readoverride.o: enums.inc dimensions.inc splitter.inc printer.inc
readoverride.o:
readparfile.o: proteins.inc dimensions.inc splitter.inc printer.inc
readparfile.o:
readpdbline.o: printer.inc
readpdbline.o:
readresidueoverride.o: enums.inc dimensions.inc splitter.inc printer.inc
readresidueoverride.o: templates.inc
readresidueoverride.o:
readresiduetemplate.o: enums.inc dimensions.inc splitter.inc constants.inc
readresiduetemplate.o: printer.inc templates.inc
readresiduetemplate.o:
readrestart.o: enums.inc proteins.inc dimensions.inc solutes.inc solvents.inc
readrestart.o: printer.inc simparams.inc flags.inc pressure.inc
readrestart.o:
readrestemplateangleline.o: enums.inc proteins.inc dimensions.inc splitter.inc
readrestemplateangleline.o: printer.inc templates.inc
readrestemplateangleline.o:
readrestemplateatomline.o: enums.inc dimensions.inc splitter.inc printer.inc
readrestemplateatomline.o: templates.inc
readrestemplateatomline.o:
readrestemplatebackboneline.o: enums.inc dimensions.inc splitter.inc
readrestemplatebackboneline.o: templates.inc
readrestemplatebackboneline.o:
readrestemplatebondline.o: enums.inc proteins.inc dimensions.inc splitter.inc
readrestemplatebondline.o: printer.inc templates.inc
readrestemplatebondline.o:
readrestemplatedihedralline.o: enums.inc proteins.inc dimensions.inc
readrestemplatedihedralline.o: splitter.inc printer.inc templates.inc
readrestemplatedihedralline.o:
readrestemplateinfoline.o: dimensions.inc templates.inc
readrestemplateinfoline.o:
readrestemplateureybradleyline.o: enums.inc dimensions.inc splitter.inc
readrestemplateureybradleyline.o: printer.inc templates.inc
readrestemplateureybradleyline.o:
readrestemplatezmatline.o: enums.inc dimensions.inc splitter.inc printer.inc
readrestemplatezmatline.o: templates.inc
readrestemplatezmatline.o:
readrotamer.o: proteins.inc dimensions.inc splitter.inc printer.inc rotamer.inc
readrotamer.o: templates.inc
readrotamer.o:
readsoltemplateangleline.o: enums.inc proteins.inc dimensions.inc splitter.inc
readsoltemplateangleline.o: printer.inc templates.inc
readsoltemplateangleline.o:
readsoltemplateatomline.o: enums.inc dimensions.inc splitter.inc printer.inc
readsoltemplateatomline.o: templates.inc
readsoltemplateatomline.o:
readsoltemplatebondline.o: enums.inc proteins.inc dimensions.inc splitter.inc
readsoltemplatebondline.o: printer.inc templates.inc
readsoltemplatebondline.o:
readsoltemplatecoordline.o: dimensions.inc printer.inc
readsoltemplatecoordline.o:
readsoltemplatedihedralline.o: enums.inc proteins.inc dimensions.inc
readsoltemplatedihedralline.o: splitter.inc printer.inc templates.inc
readsoltemplatedihedralline.o:
readsoltemplateinfoline.o: dimensions.inc templates.inc
readsoltemplateinfoline.o:
readsoltemplateureybradleyline.o: enums.inc dimensions.inc splitter.inc
readsoltemplateureybradleyline.o: printer.inc templates.inc
readsoltemplateureybradleyline.o:
readsoltemplatevariableline.o: enums.inc proteins.inc dimensions.inc
readsoltemplatevariableline.o: splitter.inc printer.inc templates.inc
readsoltemplatevariableline.o:
readsolutetemplate.o: enums.inc dimensions.inc splitter.inc printer.inc
readsolutetemplate.o: templates.inc
readsolutetemplate.o:
readsolventtemplate.o: enums.inc dimensions.inc splitter.inc printer.inc
readsolventtemplate.o: templates.inc
readsolventtemplate.o:
readsurface.o: dimensions.inc splitter.inc constants.inc printer.inc surface.inc
readsurface.o:
readsvntemplateatomline.o: enums.inc dimensions.inc splitter.inc printer.inc
readsvntemplateatomline.o: templates.inc
readsvntemplateatomline.o:
readsvntemplatecoordline.o: dimensions.inc printer.inc
readsvntemplatecoordline.o:
readsvntemplateinfoline.o: dimensions.inc templates.inc
readsvntemplateinfoline.o:
readtmpl.o: dimensions.inc splitter.inc printer.inc
readtmpl.o:
readub.o: enums.inc dimensions.inc splitter.inc parameters.inc constants.inc
readub.o: printer.inc
readub.o:
recalculategeometries.o: dimensions.inc solutes.inc
recalculategeometries.o:
recalculateparameters.o: simparams.inc dimensions.inc parameters.inc
recalculateparameters.o:
removeexcesscapsolvent.o: dimensions.inc solvents.inc
removeexcesscapsolvent.o:
resetgbenergy.o: enums.inc flags.inc dimensions.inc printer.inc
resetgbenergy.o:
residuebornintrapda.o: enums.inc flags.inc proteins.inc dimensions.inc gb.inc
residuebornintrapda.o: parameters.inc
residuebornintrapda.o:
residuedetailedintraenergy.o: nbenergy2.inc proteins.inc enums.inc
residuedetailedintraenergy.o: dimensions.inc parameters.inc constants.inc
residuedetailedintraenergy.o: nbpropro.inc printer.inc simparams.inc flags.inc
residuedetailedintraenergy.o: nbenergyvar2.inc quickpot.inc nbinitial.inc
residuedetailedintraenergy.o:
residueenergy.o: enums.inc proteins.inc dimensions.inc solutes.inc solvents.inc
residueenergy.o: constants.inc printer.inc simparams.inc energies.inc
residueenergy.o:
residueintraenergy.o: proteins.inc enums.inc dimensions.inc constants.inc
residueintraenergy.o:
residueintragbenergy.o: enums.inc proteins.inc dimensions.inc gb.inc
residueintragbenergy.o: parameters.inc constants.inc printer.inc simparams.inc
residueintragbenergy.o: flags.inc
residueintragbenergy.o:
residueintranbenergy.o: nbenergy2.inc enums.inc proteins.inc dimensions.inc
residueintranbenergy.o: parameters.inc constants.inc nbpropro.inc printer.inc
residueintranbenergy.o: nbenergyvar2.inc simparams.inc flags.inc quickpot.inc
residueintranbenergy.o: nbinitial.inc
residueintranbenergy.o:
residuemove.o: enums.inc proteins.inc dimensions.inc gb.inc constants.inc
residuemove.o: printer.inc moves.inc flags.inc energies.inc quickpot.inc
residuemove.o:
residueproteinenergy.o: enums.inc flags.inc simparams.inc proteins.inc
residueproteinenergy.o: dimensions.inc constants.inc
residueproteinenergy.o:
residueresidueborninterpda.o: enums.inc flags.inc proteins.inc dimensions.inc
residueresidueborninterpda.o: gb.inc parameters.inc printer.inc
residueresidueborninterpda.o:
residueresidueenergy.o: nbenergyvar2.inc flags.inc proteins.inc simparams.inc
residueresidueenergy.o: enums.inc dimensions.inc constants.inc
residueresidueenergy.o:
residueresidueintergbenergy.o: enums.inc proteins.inc dimensions.inc gb.inc
residueresidueintergbenergy.o: parameters.inc constants.inc printer.inc
residueresidueintergbenergy.o: simparams.inc flags.inc
residueresidueintergbenergy.o:
residuesinglesolventenergy.o: nbenergy2.inc enums.inc proteins.inc
residuesinglesolventenergy.o: dimensions.inc parameters.inc solvents.inc
residuesinglesolventenergy.o: constants.inc printer.inc nbenergyvar2.inc
residuesinglesolventenergy.o: simparams.inc flags.inc quickpot.inc nbprosvn.inc
residuesinglesolventenergy.o: nbinitial.inc
residuesinglesolventenergy.o:
residuesolutebornpda.o: enums.inc proteins.inc dimensions.inc solutes.inc gb.inc
residuesolutebornpda.o: parameters.inc constants.inc simparams.inc flags.inc
residuesolutebornpda.o:
residuesoluteenergy.o: nbenergy2.inc enums.inc proteins.inc dimensions.inc
residuesoluteenergy.o: solutes.inc nbprosol.inc parameters.inc constants.inc
residuesoluteenergy.o: printer.inc nbenergyvar2.inc simparams.inc flags.inc
residuesoluteenergy.o: quickpot.inc nbinitial.inc
residuesoluteenergy.o:
residuesolutegbenergy.o: enums.inc proteins.inc dimensions.inc solutes.inc
residuesolutegbenergy.o: gb.inc parameters.inc constants.inc printer.inc
residuesolutegbenergy.o: flags.inc simparams.inc
residuesolutegbenergy.o:
residuesolventbornpda.o: enums.inc flags.inc proteins.inc dimensions.inc gb.inc
residuesolventbornpda.o: parameters.inc solvents.inc
residuesolventbornpda.o:
residuesolventenergy.o: proteins.inc dimensions.inc constants.inc solvents.inc
residuesolventenergy.o:
residuesolventgbenergy.o: enums.inc proteins.inc dimensions.inc gb.inc
residuesolventgbenergy.o: parameters.inc solvents.inc constants.inc printer.inc
residuesolventgbenergy.o: flags.inc simparams.inc
residuesolventgbenergy.o:
residuesurface.o: enums.inc proteins.inc dimensions.inc solutes.inc
residuesurface.o: parameters.inc solvents.inc constants.inc surface.inc
residuesurface.o: printer.inc simparams.inc flags.inc
residuesurface.o:
retienergy.o: enums.inc simparams.inc dimensions.inc gb.inc constants.inc
retienergy.o: printer.inc
retienergy.o:
reversevolumemove.o: enums.inc proteins.inc dimensions.inc solutes.inc
reversevolumemove.o: solvents.inc constants.inc flags.inc simparams.inc
reversevolumemove.o: pressure.inc
reversevolumemove.o:
revertborndata.o: proteins.inc enums.inc dimensions.inc solutes.inc gb.inc
revertborndata.o: solvents.inc printer.inc flags.inc energies.inc quickpot.inc
revertborndata.o:
revertcoordinates.o: flags.inc proteins.inc enums.inc energies.inc
revertcoordinates.o: dimensions.inc quickpot.inc solutes.inc solvents.inc
revertcoordinates.o: printer.inc
revertcoordinates.o:
revertenergies.o: proteins.inc enums.inc dimensions.inc solutes.inc
revertenergies.o: constants.inc solvents.inc printer.inc flags.inc energies.inc
revertenergies.o: extraenergies.inc gb.inc surface.inc quickpot.inc
revertenergies.o:
rotamermove.o: enums.inc proteins.inc dimensions.inc constants.inc printer.inc
rotamermove.o: moves.inc flags.inc energies.inc
rotamermove.o:
rotateacid.o: proteins.inc dimensions.inc
rotateacid.o:
rotatesolute.o: dimensions.inc solutes.inc
rotatesolute.o:
rotatesolvent.o: dimensions.inc constants.inc solvents.inc
rotatesolvent.o:
saveresiduegb.o: proteins.inc dimensions.inc gb.inc printer.inc
saveresiduegb.o:
savesolutegb.o: dimensions.inc solutes.inc gb.inc printer.inc
savesolutegb.o:
savesolventgb.o: dimensions.inc gb.inc solvents.inc printer.inc
savesolventgb.o:
scanforvariablegroups.o: flags.inc simparams.inc dimensions.inc constants.inc
scanforvariablegroups.o: printer.inc
scanforvariablegroups.o:
sepresresenergy.o: nbenergy2.inc enums.inc proteins.inc dimensions.inc
sepresresenergy.o: parameters.inc constants.inc nbpropro.inc printer.inc
sepresresenergy.o: nbenergyvar2.inc flags.inc simparams.inc quickpot.inc
sepresresenergy.o: nbinitial.inc
sepresresenergy.o:
setallmoved.o: dimensions.inc movelist.inc
setallmoved.o:
setbackwardsstructure.o: simparams.inc dimensions.inc solutes.inc
setbackwardsstructure.o:
setdefaultoptions.o: enums.inc flags.inc simparams.inc dimensions.inc
setdefaultoptions.o: pressure.inc solvents.inc constants.inc printer.inc
setdefaultoptions.o:
setdirectangle.o: simparams.inc proteins.inc dimensions.inc constants.inc
setdirectangle.o: parameters.inc printer.inc
setdirectangle.o:
setdirectbond.o: simparams.inc proteins.inc dimensions.inc constants.inc
setdirectbond.o: parameters.inc printer.inc
setdirectbond.o:
setdirectdihedral.o: proteins.inc dimensions.inc parameters.inc printer.inc
setdirectdihedral.o:
setdirectureybradley.o: simparams.inc dimensions.inc constants.inc
setdirectureybradley.o: parameters.inc printer.inc
setdirectureybradley.o:
setequilibriumgeometry.o: proteins.inc enums.inc dimensions.inc solutes.inc
setequilibriumgeometry.o: parameters.inc printer.inc
setequilibriumgeometry.o:
setforwardsstructure.o: simparams.inc dimensions.inc solutes.inc
setforwardsstructure.o:
setinitialtitrationstate.o: enums.inc proteins.inc flags.inc simparams.inc
setinitialtitrationstate.o: dimensions.inc printer.inc
setinitialtitrationstate.o:
setkey.o: keyset.inc printer.inc
setkey.o:
setlength.o: dimensions.inc constants.inc
setlength.o:
setproteingeometry.o: enums.inc proteins.inc dimensions.inc printer.inc
setproteingeometry.o:
setresiduemoved.o: proteins.inc dimensions.inc printer.inc movelist.inc
setresiduemoved.o:
setsolutegeometry.o: enums.inc dimensions.inc solutes.inc printer.inc
setsolutegeometry.o:
setsolutemoved.o: dimensions.inc solutes.inc printer.inc movelist.inc
setsolutemoved.o:
setsolventmoved.o: dimensions.inc solvents.inc printer.inc movelist.inc
setsolventmoved.o:
setstream.o: dimensions.inc printer.inc
setstream.o:
setzmatrixgeometry.o: dimensions.inc printer.inc
setzmatrixgeometry.o:
shutdown.o: version.inc printer.inc
shutdown.o:
sigma.o: dimensions.inc constants.inc parameters.inc
sigma.o:
simulate.o: flags.inc enums.inc dimensions.inc pressure.inc printer.inc
simulate.o:
simulation.o: enums.inc proteins.inc dimensions.inc solutes.inc solvents.inc
simulation.o: constants.inc surface.inc rotamer.inc printer.inc averages.inc
simulation.o: simparams.inc flags.inc splitter.inc
simulation.o:
singlepointenergy.o: enums.inc dimensions.inc printer.inc
singlepointenergy.o:
solutebornintrapda.o: enums.inc flags.inc dimensions.inc solutes.inc gb.inc
solutebornintrapda.o: parameters.inc
solutebornintrapda.o:
solutedetailedintraenergy.o: nbenergy2.inc enums.inc proteins.inc dimensions.inc
solutedetailedintraenergy.o: solutes.inc parameters.inc constants.inc
solutedetailedintraenergy.o: printer.inc nbenergyvar2.inc flags.inc
solutedetailedintraenergy.o: simparams.inc nbsolsol.inc nbinitial.inc
solutedetailedintraenergy.o:
soluteenergy.o: enums.inc proteins.inc dimensions.inc solutes.inc solvents.inc
soluteenergy.o: constants.inc printer.inc flags.inc simparams.inc energies.inc
soluteenergy.o:
soluteintraenergy.o: nbenergy2.inc enums.inc dimensions.inc solutes.inc
soluteintraenergy.o: parameters.inc constants.inc printer.inc nbenergyvar2.inc
soluteintraenergy.o: simparams.inc flags.inc quickpot.inc nbsolsol.inc
soluteintraenergy.o: nbinitial.inc
soluteintraenergy.o:
soluteintragbenergy.o: enums.inc dimensions.inc solutes.inc gb.inc
soluteintragbenergy.o: parameters.inc constants.inc printer.inc simparams.inc
soluteintragbenergy.o: flags.inc
soluteintragbenergy.o:
solutemove.o: enums.inc dimensions.inc solutes.inc solvents.inc gb.inc
solutemove.o: printer.inc moves.inc flags.inc energies.inc quickpot.inc
solutemove.o:
solutesinglesolventenergy.o: nbenergy2.inc enums.inc dimensions.inc solutes.inc
solutesinglesolventenergy.o: parameters.inc constants.inc solvents.inc
solutesinglesolventenergy.o: printer.inc nbenergyvar2.inc simparams.inc
solutesinglesolventenergy.o: flags.inc quickpot.inc nbsolsvn.inc nbinitial.inc
solutesinglesolventenergy.o:
solutesoluteborninterpda.o: enums.inc dimensions.inc solutes.inc gb.inc
solutesoluteborninterpda.o: parameters.inc constants.inc simparams.inc flags.inc
solutesoluteborninterpda.o:
solutesoluteenergy.o: nbenergy2.inc enums.inc dimensions.inc solutes.inc
solutesoluteenergy.o: parameters.inc constants.inc printer.inc nbenergyvar2.inc
solutesoluteenergy.o: flags.inc simparams.inc quickpot.inc nbsolsol.inc
solutesoluteenergy.o: nbinitial.inc
solutesoluteenergy.o:
solutesoluteintergbenergy.o: enums.inc dimensions.inc solutes.inc gb.inc
solutesoluteintergbenergy.o: parameters.inc constants.inc printer.inc
solutesoluteintergbenergy.o: simparams.inc flags.inc
solutesoluteintergbenergy.o:
solutesolventbornpda.o: enums.inc dimensions.inc solutes.inc gb.inc
solutesolventbornpda.o: parameters.inc solvents.inc constants.inc simparams.inc
solutesolventbornpda.o: flags.inc
solutesolventbornpda.o:
solutesolventenergy.o: flags.inc dimensions.inc solutes.inc constants.inc
solutesolventenergy.o: solvents.inc printer.inc
solutesolventenergy.o:
solutesolventgbenergy.o: enums.inc dimensions.inc solutes.inc gb.inc
solutesolventgbenergy.o: parameters.inc solvents.inc constants.inc printer.inc
solutesolventgbenergy.o: flags.inc simparams.inc
solutesolventgbenergy.o:
solutesurface.o: enums.inc proteins.inc dimensions.inc solutes.inc
solutesurface.o: parameters.inc solvents.inc constants.inc surface.inc
solutesurface.o: printer.inc simparams.inc flags.inc
solutesurface.o:
solventbornintrapda.o: enums.inc flags.inc dimensions.inc gb.inc parameters.inc
solventbornintrapda.o: solvents.inc
solventbornintrapda.o:
solventcapenergy.o: enums.inc dimensions.inc constants.inc solvents.inc
solventcapenergy.o: printer.inc flags.inc simparams.inc extraenergies.inc gb.inc
solventcapenergy.o: surface.inc movelist.inc
solventcapenergy.o:
solventenergy.o: enums.inc proteins.inc dimensions.inc solutes.inc solvents.inc
solventenergy.o: constants.inc printer.inc flags.inc simparams.inc energies.inc
solventenergy.o:
solventintragbenergy.o: enums.inc dimensions.inc gb.inc parameters.inc
solventintragbenergy.o: solvents.inc constants.inc printer.inc simparams.inc
solventintragbenergy.o: flags.inc
solventintragbenergy.o:
solventmove.o: enums.inc dimensions.inc gb.inc solvents.inc constants.inc
solventmove.o: printer.inc moves.inc flags.inc energies.inc metropolis.inc
solventmove.o: quickpot.inc
solventmove.o:
solventsolventborninterpda.o: enums.inc flags.inc dimensions.inc gb.inc
solventsolventborninterpda.o: parameters.inc solvents.inc
solventsolventborninterpda.o:
solventsolventbulkenergy.o: flags.inc dimensions.inc constants.inc solvents.inc
solventsolventbulkenergy.o:
solventsolventenergy.o: flags.inc dimensions.inc constants.inc solvents.inc
solventsolventenergy.o:
solventsolventintergbenergy.o: enums.inc dimensions.inc gb.inc parameters.inc
solventsolventintergbenergy.o: solvents.inc constants.inc printer.inc
solventsolventintergbenergy.o: simparams.inc flags.inc
solventsolventintergbenergy.o:
solventsolventsingleenergy.o: nbenergy2.inc enums.inc dimensions.inc
solventsolventsingleenergy.o: parameters.inc constants.inc solvents.inc
solventsolventsingleenergy.o: printer.inc nbenergyvar2.inc simparams.inc
solventsolventsingleenergy.o: flags.inc quickpot.inc nbsvnsvn.inc nbinitial.inc
solventsolventsingleenergy.o:
solventsurface.o: enums.inc proteins.inc dimensions.inc solutes.inc
solventsurface.o: parameters.inc solvents.inc constants.inc surface.inc
solventsurface.o: printer.inc simparams.inc flags.inc
solventsurface.o:
specifiedangle.o: proteins.inc enums.inc dimensions.inc solutes.inc
specifiedangle.o:
specifiedbond.o: proteins.inc enums.inc dimensions.inc solutes.inc
specifiedbond.o:
split.o: dimensions.inc splitter.inc
split.o:
split2.o: dimensions.inc splitter.inc printer.inc
split2.o:
splitDDDSAsimulate.o: enums.inc dimensions.inc gb.inc surface.inc printer.inc
splitDDDSAsimulate.o: flags.inc energies.inc quickpot.inc pressure.inc
splitDDDSAsimulate.o:
splitGBSAsimulate.o: enums.inc dimensions.inc printer.inc flags.inc
splitGBSAsimulate.o: simparams.inc energies.inc extraenergies.inc gb.inc
splitGBSAsimulate.o: surface.inc quickpot.inc pressure.inc
splitGBSAsimulate.o:
splitchar.o: splitter.inc
splitchar.o:
storeborndata.o: proteins.inc enums.inc dimensions.inc solutes.inc gb.inc
storeborndata.o: solvents.inc printer.inc flags.inc energies.inc quickpot.inc
storeborndata.o:
storecoordinates.o: flags.inc proteins.inc enums.inc energies.inc dimensions.inc
storecoordinates.o: quickpot.inc solutes.inc solvents.inc printer.inc
storecoordinates.o:
storeenergies.o: proteins.inc enums.inc dimensions.inc solutes.inc constants.inc
storeenergies.o: solvents.inc printer.inc flags.inc energies.inc
storeenergies.o: extraenergies.inc gb.inc surface.inc quickpot.inc
storeenergies.o:
sumcharge.o: proteins.inc dimensions.inc constants.inc printer.inc
sumcharge.o:
sumresiduecharge.o: proteins.inc dimensions.inc constants.inc
sumresiduecharge.o:
sumsolutecharge.o: enums.inc dimensions.inc solutes.inc constants.inc
sumsolutecharge.o: printer.inc
sumsolutecharge.o:
surfradius.o: dimensions.inc constants.inc parameters.inc printer.inc
surfradius.o: surface.inc
surfradius.o:
surftension.o: dimensions.inc constants.inc parameters.inc surface.inc
surftension.o:
switchbornradii.o: proteins.inc enums.inc dimensions.inc solutes.inc gb.inc
switchbornradii.o: solvents.inc printer.inc flags.inc energies.inc quickpot.inc
switchbornradii.o:
synchfepsolutes.o: enums.inc simparams.inc dimensions.inc solutes.inc
synchfepsolutes.o: printer.inc
synchfepsolutes.o:
syncproteinzmatrix.o: proteins.inc enums.inc dimensions.inc printer.inc
syncproteinzmatrix.o:
syncsolutedummies.o: enums.inc dimensions.inc solutes.inc constants.inc
syncsolutedummies.o: printer.inc
syncsolutedummies.o:
takeacid.o: proteins.inc dimensions.inc
takeacid.o:
takebase.o: proteins.inc dimensions.inc
takebase.o:
testenergies.o: flags.inc enums.inc energies.inc dimensions.inc constants.inc
testenergies.o: printer.inc
testenergies.o:
titrationmove.o: enums.inc flags.inc proteins.inc energies.inc dimensions.inc
titrationmove.o: printer.inc moves.inc
titrationmove.o:
totalborn.o: enums.inc proteins.inc dimensions.inc solutes.inc gb.inc
totalborn.o: solvents.inc constants.inc printer.inc flags.inc simparams.inc
totalborn.o:
totalenergy.o: enums.inc proteins.inc dimensions.inc solutes.inc constants.inc
totalenergy.o: solvents.inc printer.inc flags.inc simparams.inc energies.inc
totalenergy.o: quickpot.inc movelist.inc
totalenergy.o:
totalsurface.o: enums.inc proteins.inc dimensions.inc solutes.inc solvents.inc
totalsurface.o: constants.inc printer.inc surface.inc flags.inc simparams.inc
totalsurface.o:
translateacid.o: proteins.inc dimensions.inc
translateacid.o:
translatesolute.o: dimensions.inc solutes.inc
translatesolute.o:
translatesolvent.o: dimensions.inc solvents.inc
translatesolvent.o:
updateangle.o: enums.inc proteins.inc dimensions.inc printer.inc
updateangle.o:
updatebond.o: enums.inc proteins.inc dimensions.inc printer.inc
updatebond.o:
updatebornradii.o: proteins.inc dimensions.inc solutes.inc gb.inc solvents.inc
updatebornradii.o: parameters.inc constants.inc
updatebornradii.o:
updatedihedral.o: enums.inc proteins.inc flags.inc dimensions.inc printer.inc
updatedihedral.o:
updateextraenergies.o: flags.inc energies.inc dimensions.inc quickpot.inc
updateextraenergies.o: extraenergies.inc gb.inc surface.inc constants.inc
updateextraenergies.o: printer.inc
updateextraenergies.o:
updatemultiplemovesenergies.o: proteins.inc energies.inc dimensions.inc
updatemultiplemovesenergies.o: solutes.inc solvents.inc printer.inc movelist.inc
updatemultiplemovesenergies.o: stacks.inc
updatemultiplemovesenergies.o:
updateresidueenergies.o: proteins.inc energies.inc dimensions.inc solutes.inc
updateresidueenergies.o: solvents.inc
updateresidueenergies.o:
updateresidueinternals.o: proteins.inc enums.inc dimensions.inc
updateresidueinternals.o:
updatesoluteenergies.o: proteins.inc energies.inc dimensions.inc solutes.inc
updatesoluteenergies.o: solvents.inc
updatesoluteenergies.o:
updatesoluteinternals.o: enums.inc dimensions.inc solutes.inc
updatesoluteinternals.o:
updatesolventenergies.o: proteins.inc energies.inc dimensions.inc solutes.inc
updatesolventenergies.o: solvents.inc
updatesolventenergies.o:
updatesplitmovelist.o: enums.inc dimensions.inc quickpot.inc printer.inc
updatesplitmovelist.o:
updatesplitmoveresiduelist.o: proteins.inc dimensions.inc quickpot.inc
updatesplitmoveresiduelist.o: printer.inc movelist.inc
updatesplitmoveresiduelist.o:
updatesplitmovesolutelist.o: dimensions.inc quickpot.inc solutes.inc printer.inc
updatesplitmovesolutelist.o: movelist.inc
updatesplitmovesolutelist.o:
updatesplitmovesolventlist.o: dimensions.inc quickpot.inc solvents.inc
updatesplitmovesolventlist.o: printer.inc movelist.inc
updatesplitmovesolventlist.o:
updatesurfaces.o: proteins.inc dimensions.inc solutes.inc constants.inc
updatesurfaces.o: solvents.inc printer.inc surface.inc flags.inc energies.inc
updatesurfaces.o:
updatetotalenergies.o: proteins.inc flags.inc energies.inc dimensions.inc
updatetotalenergies.o: solutes.inc solvents.inc
updatetotalenergies.o:
updateureybradley.o: enums.inc proteins.inc dimensions.inc printer.inc
updateureybradley.o:
ureybradleyenergy.o: enums.inc proteins.inc dimensions.inc solutes.inc
ureybradleyenergy.o: parameters.inc constants.inc printer.inc flags.inc
ureybradleyenergy.o:
veclength2.o: dimensions.inc constants.inc
veclength2.o:
vecnormalise.o: dimensions.inc constants.inc printer.inc
vecnormalise.o:
volumemove.o: enums.inc dimensions.inc solvents.inc printer.inc moves.inc
volumemove.o: flags.inc metropolis.inc energies.inc pressure.inc
volumemove.o:
writerestart.o: enums.inc proteins.inc dimensions.inc solutes.inc solvents.inc
writerestart.o: printer.inc simparams.inc flags.inc
writerestart.o:
zeroallpdacont.o: proteins.inc dimensions.inc solutes.inc gb.inc solvents.inc
zeroallpdacont.o: parameters.inc constants.inc
zeroallpdacont.o:
zeropdacont.o: proteins.inc dimensions.inc solutes.inc gb.inc solvents.inc
zeropdacont.o: parameters.inc constants.inc
zeropdacont.o:

