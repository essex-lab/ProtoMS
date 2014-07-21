      subroutine assigngcsoluteTemplate(igcs)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'enums.inc'
c####################################################
c
c     This routine assigns a template to gc-solute
c     igcs
c
c     MSB 1st April 09 
c
c####################################################

      integer igcs,i,nats,itmpl,iat,ipar1,ipar2
      integer findgcsoluteTemplateAtom,addCLJ,getCLJ
      integer findgcsoluteTemplate      

      if (igcs.le.0 .or. igcs.gt.NGCsolutes) return
      
      nats = NGCSAtoms(igcs)
      
      if (nats.le.0) return

      itmpl = findGCsoluteTemplate(GCSName(igcs))
      if (itmpl.le.0) then
        write(printstring,*) "Could not find a template for gc-solute ",GCSName(igcs)
        call printLine(INFO,printstring)
        call closeProgram(-1,printstring)
      endif
        
c     set the delta values for rotations and translations
      GCSADel(igcs) = GCSTemplateADel(itmpl)
      GCSRDel(igcs) = GCSTemplateRDel(itmpl)
        
c     check that there are the same numbers of atoms in the template
c     as there are in the solute
      if (NGCSAtoms(igcs).ne.NGCSTemplateAtms(itmpl)) then
        write(printstring,*) "There are not the same numbers of atoms in",
     .       " the pdb and template of solute ",GCSName(igcs)," (",
     .        NGCSAtoms(igcs),"-",NGCSTemplateAtms(itmpl),")"
        call closeProgram(-1,printstring)
      endif
      
c     assign the clj parameters
      do i=1,nats
        iat = findGCsoluteTemplateAtom(GCSAtomName(igcs,i),itmpl)
        if (iat.eq.0) then
          write(printstring,*) "Could not find atom ",GCSAtomName(igcs,i),
     .              " in the template for solute ",GCSName(igcs)
          call closeProgram(-1,printstring)
        endif
c       add the clj types for this solvent atom and assign the parameter          
        ipar1 = addCLJ(GCSTemplateAtmCLJ(itmpl,iat,1))
        ipar2 = addCLJ(GCSTemplateAtmCLJ(itmpl,iat,2))
        ZGCSPar(igcs,i) = getCLJ(ipar1,ipar2)
c       PRINT*, "Params!" ,ipar1,ipar2,ZGCSPar(igcs,i)
       enddo

c      temp hack
      return
      end
