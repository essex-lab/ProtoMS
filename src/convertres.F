      logical function convertRes(iprot,ires,iold,inew)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'enums.inc'
c#############################################################
c
c     This function converts ires of iprot, from residue
c     iold to residue inew - it is used to protonate and
c     deprotonate residues!
c
c     (C) Christopher Woods     1 February 2002
c
c##############################################################

      integer iprot,ires,i,izat
      character*4 iold,inew
      integer ipar0,ipar1,addCLJ,getCLJ
      integer itmp1,itmp2,ntmp1,ntmp2,itmat
      integer findTemplate,findTemplateAtom,findAtom
      
      convertRes = .false.
      
      if (iold.ne.ResNam(iprot,ires)) then
        write(printstring,*) "Big problem - the residue is not what you think it is!"
        call closeProgram(-1,printstring)
      endif
      
      itmp1 = findTemplate(iold)
      itmp2 = findTemplate(inew)
      if (itmp1.le.0 .or. itmp2.le.0) then
        write(printstring,*) "Problem converting residues from ",iold,inew,
     .                       " - could not find a template ",itmp1,itmp2
        call closeProgram(-1,printstring)
      endif
            
      ntmp1 = NTemplateAtms(itmp1)
      ntmp2 = NTemplateAtms(itmp2)
      
      if (ntmp1.le.0 .or. ntmp2.le.0 .or. ntmp1.ne.ntmp2) then
        write(printstring,*) "Problem converting residues from ",iold,inew,
     .             " - number of atoms in the templates is wrong! ",ntmp1,ntmp2
        call closeProgram(-1,printstring)
      endif

c     change the residue name over
      ResNam(iprot,ires) = inew
      
c     loop over all of the atoms, and change the clj parameters...
      do i=1,ntmp1
        izat = findAtom(iprot,ires,TemplateAtmName(itmp1,i,1))
        itmat = findTemplateAtom(itmp2,TemplateAtmName(itmp1,i,1))
        if (itmat.eq.0) then
          write(printstring,*) "Converting between two non-homogenous amino-acids - ",iold," ",inew
          call closeProgram(-1,printstring)
        endif
c       change the clj parameters over...
        ipar0 = addCLJ(TemplateAtmCLJ(itmp2,itmat,1))
        ipar1 = addCLJ(TemplateAtmCLJ(itmp2,itmat,2))
        ZPar(iprot,izat) = getCLJ(ipar0,ipar1)
      enddo

c     now loop over all of the bonds/angles/dihedrals and get thier parameters
c     reevaluated in light of the new atom types...
      if (NResBonds(iprot,ires).gt.0) then
        do i=1,NResBonds(iprot,ires)
          call updateBond(iprot,ires,i)
        enddo
      endif
      if (NResAngles(iprot,ires).gt.0) then
        do i=1,NResAngles(iprot,ires)
          call updateAngle(iprot,ires,i)
        enddo
      endif
      if (NResDihedrals(iprot,ires).gt.0) then
        do i=1,NResDihedrals(iprot,ires)
          call updateDihedral(iprot,ires,i)
        enddo
      endif
      if (NResUreyBradleys(iprot,ires).gt.0) then
        do i=1,NResUreyBradleys(iprot,ires)
          call updateUreyBradley(iprot,ires,i)
        enddo
      endif

      convertRes = .true.
      return
      end
