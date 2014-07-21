      
      double precision function residueDetailedIntraEnergy(iprot,ires)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'parameters.inc'
      include 'simparams.inc'
      include 'flags.inc'
      include 'printer.inc'
      include 'quickpot.inc'
c########################################################################
c
c     This is a version of residueIntraEnergy that prints out the 
c     values of all of the energy components (like soluteDetailedIntraEnergy)
c
c     (C) Christopher Woods, November 2004
c
c########################################################################

      integer iprot,ires
      double precision bndnrg,angnrg,dihnrg,ubynrg,nbnrg
      double precision bndnrgf,angnrgf,dihnrgf,ubynrgf,nbcnrgf,nbljnrgf
      double precision bndnrgb,angnrgb,dihnrgb,ubynrgb,nbcnrgb,nbljnrgb
      double precision ibndf,iangf,idihf,iubyf,ibndb,iangb,idihb,iubyb
      double precision icnrgf,icnrgb,iljnrgf,iljnrgb,imindist,inrg
      integer ic,inres
      double precision nbcnrg,nbljnrg,icnrg,iljnrg
      integer i,inbnds,inangs,indihs,inubys,ibnd,iang,idih,iuby
      double precision bondEnergy,angleEnergy,dihedralEnergy,residueResidueEnergy
      double precision ureyBradleyEnergy
      double precision cnrg,cnrgf,cnrgb,ljnrg,ljnrgf,ljnrgb
      double precision rad2Deg
      character*16 getAtomName,inam1,inam2,inam3,inam4
      include 'nbenergyvar2.inc'
      double precision charge,sigma,epslon      

c     Because there is no softcore btw residues
      double precision scl,sclf,sclb
      scl = 1.0d0
      sclf = 1.0d0
      sclb = 1.0d0
      residueDetailedIntraEnergy = ZERO
      softcore = .false.
      ddd = .false.
      if (iprot.le.0 .or. iprot.gt.NProteins .or. ires.le.0) return
      if (ires.gt.NResidues(iprot)) return
      
c     run over all bonds in the residue
      inbnds = NResBonds(iprot,ires)
      bndnrg = ZERO
      bndnrgf = ZERO
      bndnrgb = ZERO
      if (inbnds.ne.0) then
        ic = 0
        do i=1,inbnds
          ibnd = ResidueBond(iprot,ires,i)
          if (Bond(iprot,ibnd,3).ne.0) then
            ic = ic + 1
            inrg = bondEnergy(PROTEIN,iprot,ibnd,ibndb,ibndf)
            bndnrg = bndnrg + inrg
            bndnrgf = bndnrgf + ibndf
            bndnrgb = bndnrgb + ibndb
          
            inam1 = getAtomName(iprot,Bond(iprot,ibnd,1))
            inam2 = getAtomName(iprot,Bond(iprot,ibnd,2))
            write(printstring,10) inam1,inam2,
     .                            BondLength(iprot,ibnd),inrg
10          format('Bond ',a,' - ',a,' Length = ',f10.5,' Energy = ',f10.5)
            call printLine(ENERGY,printstring)
          endif
        enddo
        
        if (ic.gt.0) then
          write(printstring,*) "Total bond energy = ",bndnrg
          call prettyLine(ENERGY,printstring)
        endif
      endif
      
c     run over all of the angles in the residue...
      inangs = NResAngles(iprot,ires)
      if (inangs.ne.0) then
        ic = 0
        do i=1,inangs
          iang = ResidueAngle(iprot,ires,i)
          if (Angle(iprot,iang,4).ne.0) then
            ic = ic + 1
            inrg = angleEnergy(PROTEIN,iprot,iang,iangb,iangf)
            angnrg = angnrg + inrg
            angnrgf = angnrgf + iangf
            angnrgb = angnrgb + iangb
            
            inam1 = getAtomName(iprot,Angle(iprot,iang,1))
            inam2 = getAtomName(iprot,Angle(iprot,iang,2))
            inam3 = getAtomName(iprot,Angle(iprot,iang,3))
            
            write(printstring,11) inam1,inam2,inam3,
     .                            rad2Deg(AngleSize(iprot,iang)),inrg
11          format('Angle ',a,' - ',a,' - ',a,' Size = ',f10.5,' Energy = ',f10.5)
            call printLine(ENERGY,printstring)
          endif
        enddo
        
        if (ic.gt.0) then
          write(printstring,*) "Total angle energy = ",angnrg
          call prettyLine(ENERGY,printstring)
        endif
      endif
      
c     run over all of the dihedrals in the residue...
      indihs = NResDihedrals(iprot,ires)
      if (indihs.ne.0) then
        ic = 0
        do i=1,indihs
          idih = ResidueDihedral(iprot,ires,i)
          if (Dihedral(iprot,idih,5).ne.0) then
            ic = ic + 1
            inrg = dihedralEnergy(PROTEIN,iprot,idih,idihb,idihf)
            dihnrg = dihnrg + inrg
            dihnrgb = dihnrgb + idihb
            dihnrgf = dihnrgf + idihf
          
            inam1 = getAtomName(iprot,Dihedral(iprot,idih,1))
            inam2 = getAtomName(iprot,Dihedral(iprot,idih,2))
            inam3 = getAtomName(iprot,Dihedral(iprot,idih,3))
            inam4 = getAtomName(iprot,Dihedral(iprot,idih,4))
            
            write(printstring,12) inam1,inam2,inam3,inam4,
     .                            rad2Deg(DihedralSize(iprot,idih)),inrg
12          format('Dihedral ',a ,' - ',a,' - ',a,' - ',a,' Size = ',f10.5,' Energy = ',f10.5)
            call printLine(ENERGY,printstring)
          endif
        enddo
        
        if (ic.gt.0) then
          write(printstring,*) "Total dihedral energy = ",dihnrg
          call prettyLine(ENERGY,printstring)
        endif
      endif
      
c     run over all of the Urey-Bradley terms in the residue...
      inubys = NResUreyBradleys(iprot,ires)
      if (inubys.ne.0) then
        ic = 0
        do i=1,inubys
          iuby = ResidueUreyBradley(iprot,ires,i)
          if (UreyBradley(iprot,iuby,4).ne.0) then
            ic = ic + 1
            inrg = ureyBradleyEnergy(PROTEIN,iprot,iuby,iubyb,iubyf)
            ubynrg = ubynrg + inrg
            ubynrgf = ubynrgf + iubyf
            ubynrgb = ubynrgb + iubyb
            
            inam1 = getAtomName(iprot,UreyBradley(iprot,iuby,1))
            inam2 = getAtomName(iprot,UreyBradley(iprot,iuby,2))
            inam3 = getAtomName(iprot,UreyBradley(iprot,iuby,3))
            
            write(printstring,13) inam1,inam2,inam3,
     .                            UreyBradleySize(iprot,iuby),inrg
13          format('Urey-Bradley ',a,' - ',a,' - ',a,' Size = ',f10.5,' Energy = ',f10.5)
            call printLine(ENERGY,printstring)
          endif
        enddo
        
        if (ic.gt.0) then
          write(printstring,*) "Total Urey-Bradley energy = ",ubynrg
        endif
      endif

c     now get the intraresidue non-bonded term
c     no periodic boundaries in intramolecular NB interactions
      boundx = .false.
      boundy = .false.
      boundz = .false.
            
c     now run over all atom pairs in this residue, and calculate the 
c     intra-residue non-bonded energy...
      funciprot1 = iprot
      funciprot2 = iprot
      cnrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrg = ZERO
      ljnrgb = ZERO
      ljnrgf = ZERO
      
c     now run over all of the nb pairs in this residue      
      funciprot1 = iprot
      funciprot2 = iprot
      boundx = .false.
      boundy = .false.
      boundz = .false.
      if (QuickMove) ddd = .false.
      if (NResNB(iprot,ires).gt.0) then
        do i=1,NResNB(iprot,ires)
          funcat1 = ResNB(iprot,ires,i,1)
          funcat2 = ResNB(iprot,ires,i,2)
        
          inam1 = getAtomName(iprot,funcat1)
          inam2 = getAtomName(iprot,funcat2)
            
          include 'nbinitial.inc'
          include 'nbpropro.inc'
          include 'nbenergy2.inc'
            
          cnrg = cnrg + funccnrg
          ljnrg = ljnrg + funcljnrg
          cnrgf = cnrgf + funccnrgf
          cnrgb = cnrgb + funccnrgb
          ljnrgf = ljnrgf + funcljnrgf
          ljnrgb = ljnrgb + funcljnrgb
          write(printstring,14) 'NB  ',inam1,inam2,
     .                          dsqrt(funcmindist),funccnrg,funcljnrg
c         if the distance is less than sqrt(6.0) A, then put a box around it
          if (funcmindist.lt.6.0) then
            call prettyLine(ENERGY,printstring)
          else          
            call printLine(ENERGY,printstring)
          endif
        enddo
      endif
      
c     now the intraresidue 1-4 pairs
      if (NResNB14(iprot,ires).gt.0) then
        do i=1,NResNB14(iprot,ires)
          funcat1 = ResNB14(iprot,ires,i,1)
          funcat2 = ResNB14(iprot,ires,i,2)
          
          inam1 = getAtomName(iprot,funcat1)
          inam2 = getAtomName(iprot,funcat2)
            
          include 'nbinitial.inc'
          include 'nbpropro.inc'
          include 'nbenergy2.inc'
            
          cnrg = cnrg + OneOverScl14C*funccnrg
          ljnrg = ljnrg + OneOverScl14LJ*funcljnrg
          cnrgf = cnrgf + OneOverScl14C*funccnrgf
          cnrgb = cnrgb + OneOverScl14C*funccnrgb
          ljnrgf = ljnrgf + OneOverScl14LJ*funcljnrgf
          ljnrgb = ljnrgb + OneOverScl14LJ*funcljnrgb
          write(printstring,14) 'NB14',inam1,inam2,
     .          dsqrt(funcmindist),OneOverScl14C*funccnrg,OneOverScl14LJ*funcljnrg,
     .          ' SCLC = ',OneOverScl14C,' SCLJ = ',OneOverScl14LJ
c         if the distance is less than sqrt(6.0) A, then put a box around it
          if (funcmindist.lt.6.0) then
            call prettyLine(ENERGY,printstring)
          else          
            call printLine(ENERGY,printstring)
          endif
        enddo
      endif
14    format(a4,1x,a," - ",a,' Distance = ',f10.5,' Coul = ',g15.5,' LJ = ',g15.5,
     .       a8,f8.3,a8,f8.3 )
              
      write(printstring,*) "Total NB energy = ",cnrg+ljnrg," coul = ",cnrg," lj = ",ljnrg
      call prettyLine(ENERGY,printstring)
      
      nbnrg = cnrg + ljnrg

c     now run over all residues, and calculate the res-res energy
      inres = NResidues(iprot)
      if (inres.gt.1) then
        do i=1,inres
           if (i.ne.ires) then
             inrg = residueResidueEnergy(iprot,ires,i,
     .                       icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
             nbnrg = nbnrg + inrg
             nbcnrg = nbcnrg + icnrg
             nbljnrg = nbljnrg + iljnrg
             nbcnrgf = nbcnrgf + icnrgf
             nbljnrgf = nbljnrgf + iljnrgf
             nbcnrgb = nbcnrgb + icnrgb
             nbljnrgb = nbljnrgb + iljnrgb
             write(printstring,*) "Intra-residue-residue ",i," energy = ",inrg
             call printLine(ENERGY,printstring)
             write(printstring,*) "Coul = ",icnrg," lj = ",iljnrg
             call printLine(ENERGY,printstring)
           endif
        enddo
      endif
        
      residueDetailedIntraEnergy = bndnrg + angnrg + dihnrg + ubynrg + nbnrg
      write(printstring,*) "Total residue intra energy = ",residueDetailedIntraEnergy
      call prettyLine(ENERGY,printstring)
      
      return 
      
      end
