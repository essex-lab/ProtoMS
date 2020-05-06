      subroutine assignSolventTemplate(isvn)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'enums.inc'
c####################################################
c
c     This routine assigns a template to solvent
c     isvn
c
c     (C) Christopher Woods  5 March 2002
c
c####################################################

      integer isvn,i,nats,itmpl,iat,ipar1,ipar2
      integer findSolventTemplateAtom,addCLJ,getCLJ
      integer findSolventTemplate      

      if (isvn.le.0 .or. isvn.gt.NSolvents) return
      
      nats = NSvnAtoms(isvn)

      if (nats.le.0) return

      itmpl = findSolventTemplate(SvnName(isvn))
      if (itmpl.le.0) then
        write(printstring,*) "Could not find a template for solvent ",SvnName(isvn)
        call printLine(INFO,printstring)
        call closeProgram(-1,printstring)
      endif

c     set the delta values for rotations and translations
      SvnADel(isvn) = SvnTemplateADel(itmpl)
      SvnRDel(isvn) = SvnTemplateRDel(itmpl)
c     check that there are the same numbers of atoms in the template
c     as there are in the solvent
      if (NSvnAtoms(isvn).ne.NSvnTemplateAtms(itmpl)) then
        write(printstring,*) "There are not the same numbers of atoms in",
     .       " the pdb and template of solvent ",SvnName(isvn)," (",
     .        NSvnAtoms(isvn),"-",NSvnTemplateAtms(itmpl),")"
        call closeProgram(-1,printstring)
      endif
      
c     assign the clj parameters
      do i=1,nats
        iat = findSolventTemplateAtom(SvnAtomName(isvn,i),itmpl)
        if (iat.eq.0) then
          write(printstring,*) "Could not find atom ",SvnAtomName(isvn,i),
     .              " in the template for solvent ",SvnName(isvn)
          call closeProgram(-1,printstring)
        endif
c       add the clj types for this solvent atom and assign the parameter          
        ipar1 = addCLJ(SvnTemplateAtmCLJ(itmpl,iat,1))
        ipar2 = addCLJ(SvnTemplateAtmCLJ(itmpl,iat,2))
        ZSvnPar(isvn,i) = getCLJ(ipar1,ipar2)
      enddo

c     Assign the optimization flag
      SvnOptimization(isvn) = SvnTemplateOptimization(itmpl)

      return
      end
