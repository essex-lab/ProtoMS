      subroutine assignResidueTemplate(iprot,ires)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'templates.inc'
      include 'constants.inc'
      include 'enums.inc'
c########################################################
c
c     This subroutine assigns a template to residue ires
c     in protein iprot, based on its residue name...
c
c     (C) Christopher Woods 
c
c########################################################

      integer iprot,ires
      integer i,j,ipos,iatm,ibnd,iang,idih
      integer ipar0,ipar1
      integer itmpl,ichn,lgth,lineLength
      integer zstrt,zend,nats
      integer findResidueTemplate,getChainTemplate
      integer findAtom,zat1,zat2
      integer getCLJ,addCLJ
      integer iovr,findOverride
      logical needbuild
      double precision bndval,angval,dihval,idist
      character*4 inam,atmnam,bndnam,angnam,dihnam
      character*16 inam1,inam2,getAtomName
      integer oripos

      inam = ResNam(iprot,ires)

c     get the position of this residue in the chain      
      if (ires.eq.1 .and. ires.eq.NResidues(iprot)) then
        ipos = SINGLE
      else if (ires.eq.1) then
        ipos = FIRST
      else if (ires.eq.NResidues(iprot)) then
        ipos = LAST
      else
        ipos = MIDDLE
      endif
c     If ProteinScoop assume no Nterm/Cterm
      oripos = ipos
      if (ProteinScoop) ipos = MIDDLE

c     find the template for this residue
      itmpl = findResidueTemplate(inam)
      if (itmpl.le.0) then
        write(printstring,10) inam,ResNum(iprot,ires),iprot
10      format("Residue ",a,"(",i3,") in protein ",i1," does not",
     .          " appear to have a valid template!")
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Check your parameter file has a matching residue template")
        return
      endif      

c     now find the chain template for this residue in this position
      ichn = getChainTemplate(itmpl,ipos)
      
      if (ichn.le.0) then
        write(printstring,20) inam,ResNum(iprot,ires),iprot
20      format("Residue ",a,"(",i3,") in protein ",i1," does not have ",
     .         "a chain template for this position")
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Check your parameter file for a matching chain template")
        return
      endif

c     get the number of atoms in the residue - this is equal to the number
c     of atoms in the residue template plus the number of atoms in the 
c     chain template plus the four bbatoms
      nats = NTemplateAtms(itmpl) + NChainAtoms(ichn) + 4
      
c     the first atom of this residue is at index 'zstrt' in the ZCoords array
      if (ires.eq.1) then
        zstrt = 0
      else
c       the first atom is one above the last atom of the previous residue
c       The residues go from BBZid(iprot,ires,1) to BBZid(iprot,ires,4)
        zstrt = BBZid(iprot,ires-1,4)
      endif
      
c     the last atom is at zstrt + nats
      zend = zstrt + nats
      if (zend.gt.MAXATOMS) then
        write(printstring,*) "Exceeded maximum number of protein atoms (",
     .                       MAXATOMS,")"
        call printLine(FATAL,printstring)
        write(printstring,*)"Increase MAXATOMS and recompile"
        call closeProgram(-1,printstring)
        return
      endif      

c     store the number of atoms in the protein
      NZAtoms(iprot) = zend
      
c     make sure that all of the bbatoms in the chain template have been loaded - these
c     are critical and *must* be supplied in the pdb file!
      do i=1,4
        if (ChainBBAtom(ichn,i).eq.'   ') then
          lgth = lineLength(ChainName(ichn))
          write(printstring,30) ChainName(ichn)(1:lgth),i
          call printLine(FATAL,printstring)
30        format("Chain template ",a," is missing required bbatom ",i1,". Check your parameter file.")
          write(printstring,*)"Your simulation is now broken as this atom is very important!"
          call closeProgram(-1,printstring)
          return
        else if (BBNam(iprot,ires,i).eq.'   ') then
          lgth = lineLength(ProteinName(iprot))
          write(printstring,31) ResNam(iprot,ires),ResNum(iprot,ires),
     .                          ProteinName(iprot)(1:lgth),ChainBBAtom(ichn,i)
31        format("Residue ",a,"-",i4," of protein ",a," is missing the required bbatom ",a)
          call printLine(FATAL,printstring)
          write(printstring,*)"Your simulation is now broken as this atom is very important!"
          call closeProgram(-1,printstring)
          return
        else if (BBNam(iprot,ires,i).ne.ChainBBAtom(ichn,i)) then
          lgth = lineLength(ProteinName(iprot))
          write(printstring,32) ResNam(iprot,ires),ResNum(iprot,ires),ProteinName(iprot),
     .                          i,BBNam(iprot,ires,i),ChainBBAtom(ichn,i)
32        format("Residue ",a,"-",i4," of protein ",a," has the wrong bbatom at position ",i1,
     .             ". It has ",a," when it should have ",a)
          call printLine(FATAL,printstring)
          write(printstring,*)"Your simulation is now broken as this atom is very important!"
          call closeProgram(-1,printstring)
          return
        else
c         copy the coordinates into the ZCoords array        
          if (i.lt.4) then
            iatm = zstrt + i
          else
            iatm = zend
          endif
          
          do j=1,3
            ZCoords(iprot,iatm,j) = BBCoord(iprot,ires,i,j)
          enddo
          
c         save the index number into the ZCoords array          
          BBZid(iprot,ires,i) = iatm
          Zid2ResId(iprot,iatm,1) = ires
          Zid2ResId(iprot,iatm,2) = -i
          
c         now get the CLJ parameter indicies and assign them
          ipar0 = addCLJ(ChainBBPar(ichn,i,1))
          ipar1 = addCLJ(ChainBBPar(ichn,i,2))
          
c         now assign the clj parameter
          ZPar(iprot,iatm) = getCLJ(ipar0,ipar1)
        endif
      enddo

c     now loop over the rest of the atoms, make sure that they exist and
c     add them to the ZCoords array
      nats = NTemplateAtms(itmpl) + NChainAtoms(ichn)
      
      do i=1,nats
c       save the index number for this atom (the first three atoms are bbatoms)
        SCZid(iprot,ires,i) = zstrt + 3 + i
        
c       increment the number of sidechain atoms
        NSCAtoms(iprot,ires) = i
        
        if (i.le.NTemplateAtms(itmpl)) then
c         get the name of the atom, and of those that it is built from
          atmnam = TemplateAtmName(itmpl,i,1)
          bndnam = TemplateAtmName(itmpl,i,2)
          angnam = TemplateAtmName(itmpl,i,3)
          dihnam = TemplateAtmName(itmpl,i,4)
c         get the parameter values
c         JM comment --> Need to extend TemplateAtmCLJ so that we can have up to 6 numbers
c         For compatibility with OPLS behavior, if there were only two different parameters
c         then I should copy those two in the 4 others
c         I should use ipos to check where I lie in the protein
          if (ipos .eq. FIRST) then
            ipar0 = TemplateAtmCLJ(itmpl,i,3)
            ipar1 = TemplateAtmCLJ(itmpl,i,4)            
          else if (ipos .eq. LAST) then
            ipar0 = TemplateAtmCLJ(itmpl,i,5)
            ipar1 = TemplateAtmCLJ(itmpl,i,6)
          else 
            ipar0 = TemplateAtmCLJ(itmpl,i,1)
            ipar1 = TemplateAtmCLJ(itmpl,i,2)            
          endif
c         get the default zmatrix info from this template
          bndval = TemplateZMatVal(itmpl,i,1)
          angval = TemplateZMatVal(itmpl,i,2)
          dihval = TemplateZMatVal(itmpl,i,3)
        else
          j = i - NTemplateAtms(itmpl)
          atmnam = ChainAtm(ichn,j,1)
          bndnam = ChainAtm(ichn,j,2)
          angnam = ChainAtm(ichn,j,3)
          dihnam = ChainAtm(ichn,j,4)
          ipar0 = ChainAtmPar(ichn,j,1)
          ipar1 = ChainAtmPar(ichn,j,2)
c         get the default zmatrix info from this template
          bndval = ChainZMatVal(ichn,j,1)
          angval = ChainZMatVal(ichn,j,2)
          dihval = ChainZMatVal(ichn,j,3)
        endif
        
c       find the entry for this template atom in the loaded protein
        iatm = findAtom(iprot,ires,atmnam)
        
        needbuild = .false.
        
        if (iatm.eq.0) then
c         this atom has not been read from the PDB file - see if we could build it...
          lgth = lineLength(ProteinName(iprot))
          write(printstring,40) ResNam(iprot,ires),ResNum(iprot,ires),
     .                          ProteinName(iprot)(1:lgth),atmnam
40        format("The PDB file for residue ",a,"-",i4," of protein ",a,
     .           " is missing the template atom ",a,". Attempting to build it...")
          call printLine(WARNING,printstring)
c         set a flag to say that this atom needs to be built
          needbuild = .true.
          
c         get the index that this atom should have          
          iatm = SCZid(iprot,ires,i)
        else if (iatm.ne.SCZid(iprot,ires,i)) then
c         this is the result of a bug!        
          lgth = lineLength(ProteinName(iprot))
          write(printstring,41) ResNam(iprot,ires),ResNum(iprot,ires),
     .                          ProteinName(iprot)(1:lgth),atmnam
41        format("Severe problem with residue ",a,"-",i4," of protein ",a,
     .           " and template atom ",a,". This is the result of a bug!!!")
          call closeProgram(-1,printstring)
        endif

c       save the backward conversion from zid to res/atm id
        Zid2ResID(iprot,iatm,1) = ires
        Zid2ResID(iprot,iatm,2) = i
        
c       the atom exists, so save the coordinates!
        do j=1,3
          ZCoords(iprot,iatm,j) = SCCoord(iprot,ires,i,j)
        enddo

c       also save the parameters of this atom
        ipar0 = addCLJ(ipar0)
        ipar1 = addCLJ(ipar1)
        ZPar(iprot,iatm) = getCLJ(ipar0,ipar1)        

c       now find the indexes of the bond/angle/dihedral atoms - 
c       first the bond atom
        ibnd = findAtom(iprot,ires,bndnam)
        if (ibnd.le.0) then
          lgth = lineLength(ProteinName(iprot))
          write(printstring,42) ResNam(iprot,ires),ResNum(iprot,ires),
     .                          ProteinName(iprot)(1:lgth),"bond",atmnam,
     .                          bndnam
          call printLine(WARNING,printstring)
        endif
        ZMatSC(iprot,ires,i,1) = ibnd
        
c       now the angle atom
        iang = findAtom(iprot,ires,angnam)
        if (iang.le.0) then
          lgth = lineLength(ProteinName(iprot))
          write(printstring,42) ResNam(iprot,ires),ResNum(iprot,ires),
     .                          ProteinName(iprot)(1:lgth),"angle",atmnam,
     .                          angnam
          call printLine(WARNING,printstring)
        endif
        ZMatSC(iprot,ires,i,2) = iang
          
c       now the dihedral atom
        idih = findAtom(iprot,ires,dihnam)
        if (idih.le.0) then
          lgth = lineLength(ProteinName(iprot))
          write(printstring,42) ResNam(iprot,ires),ResNum(iprot,ires),
     .                          ProteinName(iprot)(1:lgth),"dihedral",atmnam,
     .                          dihnam
          call printLine(WARNING,printstring)
        endif        
        ZMatSC(iprot,ires,i,3) = idih

42      format("Problem with connectivity of residue ",a,"-",i4," in protein ",a,
     .         ". Cannot find the ",a," atom for atom ",a," (it should be ",a,")")

c       make sure that the connectivity is correct, e.g. iatm > ibnd > iang > idih,
c       otherwise we can't build the sidechains correctly! (except for BBAtom(4) as
c       it is the last atom of the residue, yet its coordinates are readily available)
        if (iatm.le.0 .or. ibnd.le.0 .or. iang.le.0 .or. idih.le.0) then
          if (needbuild) then
            call printLine(FATAL,"Cannot build atom as one of the required atoms has not been loaded.")
            write(printstring,*)"Your simulation is broken as you have a missing atom!"
            call closeProgram(-1,printstring)
          endif
        else if (iatm.eq.ibnd .or. iatm.eq.iang .or. iatm.eq.idih .or. ibnd.eq.iang
     .           .or. ibnd.eq.idih .or. iang.eq.idih) then
          write(printstring,44) atmnam,iatm,bndnam,ibnd,angnam,iang,dihnam,idih
44        format("Connectivity problem in template! Atom ",a,"(",i5,") should not be the same as bond atom ",
     .           a,"(",i5,") not the same as angle atom ",a,"(",i5,") and ",
     .           "not the same as dihedral atom ",a,"(",i5,")")
          call printLine(WARNING,printstring)
          call printLine(WARNING,"This is a serious problem that will break any simulation!")
          call printLine(WARNING,"Check that your residue and chain templates are correct.")
          if (needbuild) then
            call printLine(FATAL,"Cannot build atom due to this template problem.")
            write(printstring,*)"Your simulation is broken as you have a missing atom!"
            call closeProgram(-1,printstring)
          endif
        else if ( (ibnd.ne.zend .and. iatm.lt.ibnd) .or.
     .            (iang.ne.zend .and. iatm.lt.iang) .or.
     .            (idih.ne.zend .and. iatm.lt.idih) ) then
          write(printstring,43) atmnam,iatm,bndnam,ibnd,angnam,iang,dihnam,idih
43        format("Connectivity problem in template! Atom ",a,"(",i5,") should be greater than bond atom ",
     .           a,"(",i5,") and greater than angle atom ",a,"(",i5,") and ",
     .           "greater than dihedral atom ",a,"(",i5,")")
          call printLine(WARNING,printstring)
          call printLine(WARNING,"This is a serious problem that will break any simulation.")
          call printLine(WARNING,"Check that your residue and chain templates are correct.")
          if (needbuild) then
            call printLine(FATAL,"Cannot build atom due to this template problem.")
            write(printstring,*)"Your simulation is broken as you have a missing atom!"
            call closeProgram(-1,printstring)
          endif
        else if (needbuild) then
c         see if we have loaded default zmatrix values for this atom - bndval will not equal
c         zero if we have
          if (bndval.eq.ZERO) then
            call printLine(FATAL,"Cannot build this atom as there is no zmatrix information available.")
            call printLine(FATAL,"See the manual to see how to add this information.")
            write(printstring,*)"Your simulation is broken as you have a missing atom!"
            call closeProgram(-1,printstring)
          else
c           we should be able to build the atom...          
            call genProteinAtom(iprot,idih,iang,ibnd,iatm,dihval,angval,bndval)

c           save the name of the atom and set the coordinates as null                  
            SCNam(iprot,ires,i) = atmnam

c           copy the coordinates into the SCCoord array          
            do j=1,3
              SCCoord(iprot,ires,i,j) = ZCoords(iprot,iatm,j)
            enddo
            call printLine(WARNING,
     .           "...the atom has been built. Please check that the atom is in the right position")
          endif
        endif
      enddo

c     now that all of the atoms have been located, see if there are any override parameters
c     that need to be assigned to this residue in this position - first see if there
c     are overrides for this residue in any position...
      iovr = findOverride(inam,ALLPOS)
      if (iovr.gt.0) call assignResidueOverrides(iprot,ires,iovr)
      
c     now do this residue in this position
      iovr = findOverride(inam,ipos)
      if (iovr.gt.0) call assignResidueOverrides(iprot,ires,iovr)
      
c     first we add bonds between residues - the bond goes from bbatom 1 of this residue
c     to bbatom 3 of the previous residue - only need to do this is we are not the first residue
      if (oripos.eq.LAST .or. oripos.eq.MIDDLE) then
        zat1 = BBZid(iprot,ires-1,3)
        zat2 = BBZid(iprot,ires,1)
c       CHECK THE DISTANCE BETWEEN ATOMS AND IF THE DISTANCE IS GREATER THAN 3 A THEN
C       ONLY ADD A DUMMY BOND, NOT A REAL BOND - THIS ALLOWS CORRECT TREATMENT OF 
C       SCOOPED PROTEINS...
        idist = dsqrt((ZCoords(iprot,zat1,1)-ZCoords(iprot,zat2,1))**2 + 
     .                (ZCoords(iprot,zat1,2)-ZCoords(iprot,zat2,2))**2 + 
     .                (ZCoords(iprot,zat1,3)-ZCoords(iprot,zat2,3))**2)
     
        if (idist.gt.THREE) then
c         add the interresidue bond as a dummy bond
          call addProteinBond(iprot,zat1,zat2,.true.,.false.,ZERO,-1,-1)
          
          inam1 = getAtomName(iprot,zat1)
          inam2 = getAtomName(iprot,zat2)
          lgth = lineLength(ProteinName(iprot))
          
          write(printstring,50) iprot,ProteinName(iprot)(1:lgth),inam1,inam2,idist
          call printLine(WARNING,printstring)
50        format("Backbone atoms of protein ",i3,"(",a,") are too far apart. The distance between atoms ",
     .           a," and ",a," is ",f8.3," A. Will now assume that this is a scoop, so will add a ",
     .           "dummy bond between these bonds. Do not move the backbone of these residues!")
        else
c         add a non-dummy, non-flexible bond with no overridden parameters        
          call addProteinBond(iprot,zat1,zat2,.false.,.false.,ZERO,-1,-1)
        endif
      endif
                                           
c     Now we add the chain specified bonds/angles/dihedrals/ubys
      if (NChainBonds(ichn).gt.0) then
        do i=1,NChainBonds(ichn)
          call addChainBond(iprot,ires,ichn,i)
        enddo
      endif
      
      if (NChainAngles(ichn).gt.0) then
        do i=1,NChainAngles(ichn)
          call addChainAngle(iprot,ires,ichn,i)
        enddo
      endif
      
      if (NChainDihedrals(ichn).gt.0) then
        do i=1,NChainDihedrals(ichn)
          call addChainDihedral(iprot,ires,ichn,i)
        enddo
      endif
      
      if (NChainUreyBradleys(ichn).gt.0) then
        do i=1,NChainUreyBradleys(ichn) 
          call addChainUreyBradley(iprot,ires,ichn,i)
        enddo
      endif
                  
c     Now add all of the template specified bonds.
c     The addbond routine will make sure that the bonds actually
c     exist in the residue. The template must contain all of the
c     bonds that exist within this residue
      if (NTemplateBonds(itmpl).gt.0) then
        do i=1,NTemplateBonds(itmpl)
          call addBond(iprot,ires,itmpl,i)        
        enddo
      endif
      
c     now do the same for the angles...
      if (NTemplateAngles(itmpl).gt.0) then
        do i=1,NTemplateAngles(itmpl)
           call addAngle(iprot,ires,itmpl,i)        
        enddo
      endif
 
c     now do the same for the dihedrals
      if (NTemplateDihedrals(itmpl).gt.0) then
        do i=1,NTemplateDihedrals(itmpl)
          call addDihedral(iprot,ires,itmpl,i)
        enddo
      endif

c     now the Urey-Bradley terms
      if (NTemplateUreyBradleys(itmpl).gt.0) then
        do i=1,NTemplateUreyBradleys(itmpl)
          call addUreyBradley(iprot,ires,itmpl,i)
        enddo
      endif
c     now assign any rotamer states from the rotamer library entry for this residue
c     (this catches proteins loaded after the rotamer library has been loaded - not
c     that that is possible in this code...)
      call assignRotamerLibrary(iprot,ires)

c     finally, assign the maximum translation and rotation amounts
c     for this residue
      ResADel(iprot,ires) = TemplateResADel(itmpl)
      ResRDel(iprot,ires) = TemplateResRDel(itmpl)      
      
      return
      
      end
