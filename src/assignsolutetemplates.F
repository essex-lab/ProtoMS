      subroutine assignSoluteTemplates(isol)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'constants.inc'
      include 'templates.inc'
      include 'solutes.inc'
      include 'simparams.inc'
      include 'enums.inc'
c########################################################
c
c     This routine assigns the molecule template to 
c     solute isol, based on its molecule name
c
c     (C) Christopher Woods   18 February 2002
c
c     Modified by Samuel Genheden, 2012-2013
c
c########################################################

      integer isol,tmol,findMoleculeTemplate,i,j,newsol
      double precision tempval,sumSoluteCharge
      double precision templam
      double precision charge,sigma,epslon
      character*2 amber
      character*4 oldres,ires
      character*4 iatm
      character*10 inam1,inam2,inam3,inam4,getSolAtomName
      integer nres 
      integer duplicateSolute
      integer cljID
      integer lgth,lineLength
      logical ok

      if (isol.le.0 .or. isol.gt.NSolutes) return
      
      lgth = lineLength(SoluteName(isol))
      write(printstring,*) "Fitting molecule templates to ",SoluteName(isol)(1:lgth)
      call printLine(INFO,printstring)
      tmol = findMoleculeTemplate(SoluteName(isol))
      if (tmol.le.0) then
        write(printstring,10) isol,SoluteName(isol)(1:lgth)
10      format("Solute ",i3," (",a,") does not appear to have a valid template!")
        call printLine(WARNING,printstring)
        call printLine(WARNING,
     .       "Ensure that there is a solute template with this name in one of your parameter files")
        call prettyLine(WARNING,"ProtoMS is continuing, but the simulation is now very broken!")
        return
      endif

c     assign the maxmimum translation/rotation of the solute
      SolADel(isol) = MolTemplateADel(tmol)
      SolRDel(isol) = MolTemplateRDel(tmol)   

c     SG 2012-2013, assign JAWS-specific parameters

c     assign the maximum theta movement of the solute
      SolTDel(isol) = MolTemplateTDel(tmol)

c     assign JAWS parameters of the solute
      SolJCorr(isol) = MolTemplateJCorr(tmol)
      SolJPmfOrd(isol) = MolTemplateJPmfOrd(tmol)
      do i=1,SolJPmfOrd(isol)
        SolJPmfCoeff(isol,i) = MolTemplateJPmfCoeff(tmol,i)
      enddo
      
c     the number of atoms in the solute must be the same as the number of atoms in the template      
      NSolAtoms(isol) = NMolTemplateAtms(tmol)

c     now that all of the coordinates have been assigned, add the solute dummy atoms
      call addSoluteDummies(isol)
      
c     loop over all of the template atoms and set the residue names and numbers of
c     the loaded atoms - also ensure that all atoms have been loaded!
      if (NMolTemplateAtms(tmol).gt.0) then
        nres = 0
        do i=1,NMolTemplateAtms(tmol)
          iatm = MolTemplateAtmName(tmol,i)
          ires = MolTemplateAtmRes(tmol,i)

c         save the residue numbers          
          if (nres.eq.0) then
            nres = 1
            oldres = ires
            SolResNam(isol,nres) = ires
c           save the residue number (this is a bit useless! It is a carry-over 
c           from when the residue number was obtained from the PDB, and is staying
c           to remind me that I need to get the residue number from the PDB)            
            SolResNum(isol,nres) = nres
            SolResidues(isol,nres,1) = i
          else if (ires.ne.oldres) then
c           save the range of atoms in each residue            
            SolResidues(isol,nres,2) = i-1
            nres = nres + 1
            if (nres.gt.MAXSOLUTERESIDUES) then
              write(printstring,*) "Exceeded maximum number of residues per solute (",
     .                           MAXSOLUTERESIDUES,")"
              call printLine(WARNING,printstring)
              call printLine(WARNING,"Increase MAXSOLUTERESIDUES and recompile")
              nres = MAXSOLUTERESIDUES
            endif
            
            SolResidues(isol,nres,1) = i
            oldres = ires
            
            SolResNam(isol,nres) = ires
            SolResNum(isol,nres) = nres
          endif
          
          SolResidue(isol,i) = nres
                    
c         the solute has been loaded such that if this atom existed in the PDB, it
c         has been loaded at this index - check this...
          if (iatm .ne. SolAtomNam(isol,i)) then
            write(printstring,*) "Atom ",iatm," in residue ",ires," was not found in the PDB. ",
     .                         "Will attempt to build this atom (please check the coordinates)"
            call printLine(WARNING,printstring)

            call printLine(WARNING,"Uh-oh! This version of the code does not yet support the building")
            call printLine(WARNING,"of missing solute atoms! This feature is planned though...")
            call printLine(WARNING,"Until then, this simulation is now very broken!")
            
            SolAtomNam(isol,i) = iatm
            SolAtomNum(isol,i) = 9999
            ZSolCoords(isol,i,1) = ZERO
            ZSolCoords(isol,i,2) = ZERO
            ZSolCoords(isol,i,3) = ZERO
          endif
        enddo
c       save the last atom number of the last residue
        SolResidues(isol,nres,2) = NMolTemplateAtms(tmol)
        NSolResidues(isol) = nres
      endif    

c     loop over the template atoms and add them (gives parameters, and
c     constructs the zmatrix!)
      if (NMolTemplateAtms(tmol).gt.0) then 
        do i=1,NMolTemplateAtms(tmol)
          call addMolAtom(isol,tmol,i)
        enddo
      endif
      
c     it is *very* important that the entire solute zmatrix is assigned!
c     thus we will check it here...
      ok = .true.
      do i=1,NSolAtoms(isol)
        do j=1,3
          if (ZSolMat(isol,i,j).eq.0) then
            write(printstring,*) "The zmatrix has not been completely assigned! ",
     .             "Atom ",i,"(",SolAtomNam(isol,i),"-",SolResidue(isol,i),
     .             ") is missing a ",j," (1=bond, 2=angle, 3=dihedral)"
            call printLine(WARNING,printstring)
            ok = .false.
          endif
        enddo
      enddo
      if (.not.ok) then
        write(printstring,*) "There have been errors setting up the zmatrix for solute ",
     .                       SoluteName(isol)(1:lgth)
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Check the pdb and template input for this solute")
        call prettyLine(WARNING,"ProtoMS will continue but the simulation is now broken!")
        return
      endif
      
c     now that we have the zmatrix, calculate all the internal coordinates
c     necessary from the pdb coordinates
      call getSoluteZmatrix(isol)
      
c     now check the charges on the solute...
      tempval = sumSoluteCharge(isol)      
      
c     loop over the template bonds and add them
      if (NMolTemplateBonds(tmol).gt.0) then
        do i=1,NMolTemplateBonds(tmol)        
          call addMolBond(isol,tmol,i)
        enddo
      endif
      
c     now the angles
      if (NMolTemplateAngles(tmol).gt.0) then
        do i=1,NMolTemplateAngles(tmol)
          call addMolAngle(isol,tmol,i)
        enddo
      endif
      
c     now the dihedrals
      if (NMolTemplateDihedrals(tmol).gt.0) then
        do i=1,NMolTemplateDihedrals(tmol)
          call addMolDihedral(isol,tmol,i)
        enddo
      endif

c     now the urey-bradley terms
      if (NMolTemplateUreyBradleys(tmol).gt.0) then
        do i=1,NMolTemplateUreyBradleys(tmol)
          call addMolUreyBradley(isol,tmol,i)
        enddo
      endif

c     now get all of the non-bonded pairs in this solute
      call getSoluteNonBonded(isol)      
      
c     assign all of the variable geometries from the template
      if (NMolTemplVar(tmol).gt.0) then
        write(printstring,*) "Building geometry duplicates for ",SoluteName(isol)
        call printLine(INFO,printstring)
c       assign the variable atoms
        call assignSoluteVariables(isol,tmol)
        if (NSolVariable(isol).le.0) then
          call closeProgram(-1,"Strange, shouldn't there be variable geometries?")
        endif
      endif

c     now clean up the solute - remove all dummy/zero bonds/angles/dihedrals
      call cleanSoluteParameters(isol)

c     now everything has been done, rebuild the solute to make sure everything is calculated
c     and initialised properly
      call buildSolute(isol)

c     now assign variable geometries, and set the forwards and backwards structures
c     to lambda = 1 and lambda = 0. This is necessary to get the non-bonded list 
c     calculation right
      if (NSolVariable(isol).gt.0) then
c       duplicate the solute
        newsol = duplicateSolute(isol)
        if (newsol.le.0) then
          write(printstring,*) "Problem duplicating solute ",SoluteName(isol)
          call closeProgram(-1,printstring)
        endif
c       set this as the forwards structure - built with lambda=1
        templam = lambdaF
        lambdaF = ONE
        call setForwardsStructure(isol,newsol,tmol)
        lambdaF = templam
c       do this again for the backwards structure
        newsol = duplicateSolute(isol)
        if (newsol.le.0) then
          write(printstring,*) "Problem duplicating solute ",SoluteName(isol)
          call closeProgram(-1,printstring)
        endif
c       set this as the backwards structure, with lambda=0
        templam = lambdaB
        lambdaB = ZERO
        call setBackwardsStructure(isol,newsol,tmol)
        lambdaB = templam
      endif
      
c     now build the correct geometry duplicates for the reference, forwards and backwards 
c     values of lambda...
      if (NSolVariable(isol).gt.0) then
        call synchFEPSolutes(isol)
c       we need to fix things so that the non-bonded list is copied to the forwards
c       and backwards solutes as well...
        do i=1,NSolNB(isol)
          do j=1,2
            SolNB(NSolForwards(isol),i,j) = SolNB(isol,i,j)
            SolNB(NSolBackwards(isol),i,j) = SolNB(isol,i,j)
          enddo
          SolNB14(NSolForwards(isol),i) = SolNB14(isol,i)
          SolNB14(NSolBackwards(isol),i) = SolNB14(isol,i)
        enddo
      endif
      
c     now print out the information about all of the assigned atoms
c     bonds/angles/dihedrals
      write(printstring,20) "ZID","ATM","RES","RES","X","Y","Z","INTPR","PAR","AM",
     .                      "CHG","SIG","EPS"
      call printLine(DETAIL,printstring)
      do i=1,NSolAtoms(isol)
        write(printstring,11) i,SolAtomNam(isol,i),SolResidue(isol,i),
     .             SolResNam(isol,SolResidue(isol,i)),ZSolCoords(isol,i,1),
     .             ZSolCoords(isol,i,2),ZSolCoords(isol,i,3),ZSolPar(isol,i),
     .             cljID(ZSolPar(isol,i)),amber(ZSolPar(isol,i)),
     .             charge(ZSolPar(isol,i)),
     .             sigma(ZSolPar(isol,i)),
     .             epslon(ZSolPar(isol,i))
        call printLine(DETAIL,printstring)
      enddo
11    format(i5,1x,a,1x,i5,1x,a,1x,3(f8.3,1x),i5,1x,i5,1x,a2,1x,3(f8.3,1x))
20    format(1x,a,2x,a,2x,1x,a,2x,a,1x,3(4x,a,4x),a,2x,a,2x,a,1x,3(3x,a,3x))

      if (NSolBonds(isol).gt.0) then
        do i=1,NSolBonds(isol)
          if (SolBond(isol,i,3).ne.0) then
            inam1 = getSolAtomName(isol,SolBond(isol,i,1))
            inam2 = getSolAtomName(isol,SolBond(isol,i,2))
            write(printstring,12) i,inam1,inam2,SolBond(isol,i,3)
            call printLine(DETAIL,printstring)
          endif
        enddo
12      format("Bond ",i5,2x,a," - ",a,1x,i5)
      endif
      
      if (NSolAngles(isol).gt.0) then
        do i=1,NSolAngles(isol)
          if (SolAngle(isol,i,4).ne.0) then
            inam1 = getSolAtomName(isol,SolAngle(isol,i,1))
            inam2 = getSolAtomName(isol,SolAngle(isol,i,2))
            inam3 = getSolAtomName(isol,SolAngle(isol,i,3))
            write(printstring,14) i,inam1,inam2,inam3,SolAngle(isol,i,4)
            call printLine(DETAIL,printstring)
          endif
        enddo
      endif
14    format("Angle ",i5,2x,a," - ",a," - ",a,1x,i5)
      
      if (NSolDihedrals(isol).gt.0) then
        do i=1,NSolDihedrals(isol)
          if (SolDihedral(isol,i,5).ne.0) then
            inam1 = getSolAtomName(isol,SolDihedral(isol,i,1))
            inam2 = getSolAtomName(isol,SolDihedral(isol,i,2))
            inam3 = getSolAtomName(isol,SolDihedral(isol,i,3))
            inam4 = getSolAtomName(isol,SolDihedral(isol,i,4))
            write(printstring,15) i,inam1,inam2,inam3,inam4,SolDihedral(isol,i,5)
            call printLine(DETAIL,printstring)
          endif
        enddo
      endif
15    format("Dihedral ",i5,2x,a,3(" - ",a),1x,i5)

      if (NSolUreyBradleys(isol).gt.0) then
        do i=1,NSolUreyBradleys(isol)
          if (SolUreyBradley(isol,i,4).ne.0) then
            inam1 = getSolAtomName(isol,SolUreyBradley(isol,i,1))
            inam2 = getSolAtomName(isol,SolUreyBradley(isol,i,2))
            inam3 = getSolAtomName(isol,SolUreyBradley(isol,i,3))
            write(printstring,16) i,inam1,inam2,inam3,SolUreyBradley(isol,i,4)
            call printLine(DETAIL,printstring)
          endif
        enddo
16      format("Urey-Bradley ",i5,2x,a," - ",a," - ",a,1x,i5)
      endif
         
      write(printstring,*) "Finished fitting template..."
      call printLine(INFO,printstring)
               
      return
      end
