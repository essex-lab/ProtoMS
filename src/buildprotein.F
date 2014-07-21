      subroutine buildProtein(iprot)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'proteins.inc'
c##################################################
c
c     This subroutine builds the entire protein,
c     iprot using coordinates from the pdb file
c     and connectivity from the templates
c
c     (C) Christopher Woods, October 2004
c
c##################################################

      integer iprot
      integer i,j,k,nres,iclj,nats,iatm
      integer cljID
      character*16 inam1,inam2,inam3,inam4,getAtomName
      character*2 amber
      double precision rad2Deg,charge,sigma,epslon
      integer inflx

      if ((iprot.gt.MAXPROTEINS).or.(iprot.le.0)) then
        write(printstring,*) "Requesting the building of an invalid protein! (",iprot,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is probably a programming bug!")
        return
      endif
      
      nres = NResidues(iprot)
      
      if (nres.le.0) then
        write(printstring,*) "Trying to build a protein that contains no residues!"
        call printLine(WARNING,printstring)
        return
      endif
      
      call printLine(INFO,"Fitting amino acids to a template...")

c     assign all of the residue templates - this also builds missing atoms            
      do i=1,nres
        call assignResidueTemplate(iprot,i)
      enddo

c     now that all of the bonds/angles/dihedrals are added, generate the 
c     residue-residue bond list...
      call proteinResResBondList(iprot)
      
c     now all of the atoms exist and all of the interatomic bonds exist
c     we can get the non-bonded list and add any missing implicit
c     angles/dihedrals
      call getNonBonded(iprot)

c     now sync the protein zmatrix - this ensures that the entire
c     data structure is internally consistent with the pdb coordinates
      call syncProteinZmatrix(iprot)
                  
c     now print the total charges on the protein..
      call sumCharge(iprot)

c     JM Mod . Initialise the number of flexible residues in this protein
c     This can be modified if a chunk latter fixsome residues
c     now work out how many flexible residues there are in this protein

      inflx = 0
      do i=1,NResidues(iprot)
        if (.not.fixedRes(iprot,i)) then 
          inflx = inflx + 1
          ResMove(iprot,inflx) = i
        endif
      enddo
      NResMove(iprot) = inflx

c     write out some statistics about the protein
      write(printstring,*) "Protein contains ",NZAtoms(iprot)," atoms, ",NBonds(iprot)," bonds, ",
     .   NAngles(iprot)," angles, ",NDihedrals(iprot)," dihedrals and ",NUreyBradleys(iprot),
     .    " Urey Bradley terms."
      call printLine(INFO,printstring)     

c     if we have 'DETAIL' output then output the entire protein zmatrix
      if (streamhandle(DETAIL).lt.0) return

      do i=1,nres
        write(printstring,*) "Residue ",ResNam(iprot,i)," (",ResNum(iprot,i),")"
        call printLine(DETAIL,printstring)
        
c       print out the backbone atoms of this residue        
        do j=1,4
          iclj = ZPar(iprot,BBZid(iprot,i,j))
          write(printstring,10) j,BBZid(iprot,i,j),BBNam(iprot,i,j),cljID(iclj),
     .         (BBCoord(iprot,i,j,k),k=1,3),charge(iclj),
     .          amber(iclj),sigma(iclj),epslon(iclj)
10        format("BBAtom ",i1,1x,i5,1x,a,1x,i5,4(1x,f8.3),1x,a,2(1x,f8.3))
          call printLine(DETAIL,printstring)
        enddo

c       now print out the sidechain atoms                
        nats = NSCAtoms(iprot,i)
        do j=1,nats
          iatm = SCZid(iprot,i,j)
          iclj = ZPar(iprot,iatm)
          
          write(printstring,11) j,iatm,SCNam(iprot,i,j),cljID(iclj),
     .          ZMatSC(iprot,i,j,1),ZValSC(iprot,i,j,1),ZMovSC(iprot,i,j,1),
     .    (ZMatSC(iprot,i,j,k),rad2Deg(ZValSC(iprot,i,j,k)),ZMovSC(iprot,i,j,k),k=2,3),
     .    charge(iclj),amber(iclj),sigma(iclj),epslon(iclj)
11        format("SCAtom ",i3,1x,i5,1x,a,1x,i5,3(1x,i5,1x,f8.3,1x,l1),1x,f8.3,1x,a,
     .                     2(1x,f8.3))
          call printLine(DETAIL,printstring)
        enddo

c       now print out all of the bonds/angles/dihedrals/ubys that connect to this
c       residue (inter-residue terms will be double-printed, but this is a 
c       small price to pay for legibility!)        
        nats = NResBonds(iprot,i)
        do j=1,nats
          iatm = ResidueBond(iprot,i,j)
          inam1 = getAtomName(iprot,Bond(iprot,iatm,1))
          inam2 = getAtomName(iprot,Bond(iprot,iatm,2))
          write(printstring,12) j,inam1,inam2,
     .                            Bond(iprot,iatm,3),BondDummy(iprot,iatm)                
12        format("Bond ",i3,1x,2(a,1x),i5,1x,l1)
          call printLine(DETAIL,printstring)
        enddo

        nats = NResAngles(iprot,i)
        do j=1,nats
          iatm = ResidueAngle(iprot,i,j)
          inam1 = getAtomNAme(iprot,Angle(iprot,iatm,1))
          inam2 = getAtomNAme(iprot,Angle(iprot,iatm,2))
          inam3 = getAtomNAme(iprot,Angle(iprot,iatm,3))
          write(printstring,13) j,inam1,inam2,inam3,
     .                            Angle(iprot,iatm,4),AngleDummy(iprot,iatm)                
13        format("Angle ",i3,1x,3(a,1x),i5,1x,l1)
          call printLine(DETAIL,printstring)
        enddo
             
        nats = NResDihedrals(iprot,i)
        do j=1,nats
          iatm = ResidueDihedral(iprot,i,j)
          inam1 = getAtomNAme(iprot,Dihedral(iprot,iatm,1))
          inam2 = getAtomNAme(iprot,Dihedral(iprot,iatm,2))
          inam3 = getAtomNAme(iprot,Dihedral(iprot,iatm,3))
          inam4 = getAtomNAme(iprot,Dihedral(iprot,iatm,4))
          write(printstring,14) j,inam1,inam2,inam3,inam4,
     .                            Dihedral(iprot,iatm,5),DihedralDummy(iprot,iatm)                
14        format("Dihedral ",i3,1x,4(a,1x),i5,1x,l1)
          call printLine(DETAIL,printstring)
        enddo
             
        nats = NResUreyBradleys(iprot,i)
        do j=1,nats
          iatm = ResidueUreyBradley(iprot,i,j)
          inam1 = getAtomNAme(iprot,UreyBradley(iprot,iatm,1))
          inam2 = getAtomNAme(iprot,UreyBradley(iprot,iatm,2))
          inam3 = getAtomNAme(iprot,UreyBradley(iprot,iatm,3))
          write(printstring,15) j,inam1,inam2,inam3,UreyBradley(iprot,iatm,4)
15        format("Urey-Bradley ",i3,1x,3(a,1x),i5)
          call printLine(DETAIL,printstring)
        enddo
        
      enddo
            
      end
