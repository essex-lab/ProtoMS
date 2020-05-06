      subroutine assignSoluteVariables(isol,tmol)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'solutes.inc'
      include 'simparams.inc'
c#########################################################
c
c     This assigns all of the variable atoms in the
c     template tmol, to solute isol
c
c     (C) Christopher Woods   14 August 2002
c
c#########################################################

      integer isol,zatm,tmol,i
      integer findSoluteAtom          
      double precision perturbForceField
      character*4 iatm,ires

c     loop through all of the variable geometries, and set
c     them in the molecule
      do i=1,NMolTemplVar(tmol)         
        iatm = MolTemplVarNam(tmol,i)
        ires = MolTemplVarRes(tmol,i)
          
        zatm = findSoluteAtom(isol,iatm,ires)        
        if (zatm.eq.0) then
          write(printstring,*) "Could not locate variable atom ",iatm,"-",ires
          call printLine(WARNING,printstring)
          return
        else if (zatm.lt.0) then
          write(printstring,*) "Cannot make dummy atmom ",iatm,"-",ires," vary geometry with lambda!"
          call printLine(WARNING,printstring)
          return
        endif                             
        NSolVariable(isol) = NSolVariable(isol) + 1
        if (NSolVariable(isol).gt.MAXSOLVARIABLES) then
          write(printstring,*) "Exceeded maximum number of variable geometries per ",
     .           "solute ",MAXSOLVARIABLES
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXSOLUTEVARIABLES and recompile"
          call closeProgram(-1,printstring)
          NSolVariable(isol) = NSolVariable(isol) - 1
          return
        endif
         
        SolVariableAtm(isol,NSolVariable(isol)) = zatm
        SolVariableTyp(isol,NSolVariable(isol)) = MolTemplVarTyp(tmol,i)
        SolVariableVal(isol,NSolVariable(isol),1) = MolTemplVarVal(tmol,i,1)
        SolVariableVal(isol,NSolVariable(isol),2) = MolTemplVarVal(tmol,i,2)
          
c       now set the value for lambda...
        ZSolVal(isol,SolVariableAtm(isol,NSolVariable(isol)),
     .                     SolVariableTyp(isol,NSolVariable(isol))) = 
     .  perturbForceField(SolVariableVal(isol,NSolVariable(isol),1),
     .                     SolVariableVal(isol,NSolVariable(isol),2),
     .                     lambda)
          
      enddo

      return
      end
