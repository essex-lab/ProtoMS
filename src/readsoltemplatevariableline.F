      subroutine readSolTemplateVariableLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c#############################################################
c
c     This routine reads in the molecular variations for
c     free energy changes associated with geometric variations
c
c##############################################################

      character*(*) line
      integer itmpl,ivar
      character*4 iatm,ires
      double precision ival0,ival1
      double precision readDouble,deg2Rad
           
c     format is 'variable iatm ires bond|angle|dihedral ival0 ival1'
      if (nwords.lt.6) then
        call printLine(WARNING,"Incorrectly formatted solute variable geometry template line?")
        call printLine(WARNING,line)
        return
      endif
      
c     increment the number of variable geometries
      NMolTemplVar(itmpl) = NMolTemplVar(itmpl) + 1
      ivar = NMolTemplVar(itmpl)
      if (ivar.gt.MAXMOLTEMPLVAR) then
        write(printstring,"(a,i8,a)")
     .        "Exceeded maximum number of template solute variations (",MAXMOLTEMPLVAR,")"
        call printLine(FATAL,printstring)
        call closeProgram(-1,"Increase MAXMOLTEMPVAR and recompile")
        NMolTemplVar(itmpl) = NMolTemplVar(itmpl) - 1
      endif
      
c     read in the atm, residue
      iatm = words(2)(1:4)
      ires = words(3)(1:4)

c     save the name and residue
      MolTemplVarNam(itmpl,ivar) = iatm
      MolTemplVarRes(itmpl,ivar) = ires

c     get the value of the geometry at lam=0 and lam=1
      ival0 = readDouble(words(5))
      ival1 = readDouble(words(6))
      
c     are we looking at a bond/angle/dihedral?
10    format(2(a,i8),2(a,f12.4))
      if (words(4).eq.'bond') then
        write(printstring,10)
     .        "Variable bond for atom ",iatm,"-",ires," starting from ",ival0," going to ",ival1
        call printLine(INFO,printstring)
        MolTemplVarTyp(itmpl,ivar) = BONDS
        MolTemplVarVal(itmpl,ivar,1) = ival0
        MolTemplVarVal(itmpl,ivar,2) = ival1        
      else if (words(4).eq.'angle') then
c       get the smallest perturbation
        call getSmallestDelta(ival0,ival1)
        write(printstring,10)
     .        "Variable angle for atom ",iatm,"-",ires," starting from ",ival0," going to ",ival1
        call printLine(INFO,printstring)
        MolTemplVarTyp(itmpl,ivar) = ANGLES
        MolTemplVarVal(itmpl,ivar,1) = deg2Rad(ival0)
        MolTemplVarVal(itmpl,ivar,2) = deg2Rad(ival1)
      else if (words(4).eq.'dihedral') then
c       we need to ensure that the perturb over the smallest delta
        call getSmallestDelta(ival0,ival1)
        write(printstring,10)
     .        "Variable dihedral for atom ",iatm,"-",ires," starting from ",ival0," going to ",ival1
        call printLine(INFO,printstring)
        MolTemplVarTyp(itmpl,ivar) = DIHEDRALS
        MolTemplVarVal(itmpl,ivar,1) = deg2Rad(ival0)
        MolTemplVarVal(itmpl,ivar,2) = deg2Rad(ival1)
      else
        call printLine(WARNING,"Cannot recognise solute variation type?")
        call printLine(WARNING,line)
        NMolTemplVar(NMolTemplates) = NMolTemplVar(NMolTemplates) - 1
        return
      endif

      end
