      double precision function residueIntraEnergy(iprot,ires,
     .                   bndnrg, angnrg, dihnrg, ubynrg, nbcnrg, nbljnrg,
     .                   bndnrgb,angnrgb,dihnrgb,ubynrgb,nbcnrgb,nbljnrgb,
     .                   bndnrgf,angnrgf,dihnrgf,ubynrgf,nbcnrgf,nbljnrgf)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'printer.inc'
c########################################################################
c
c     This function will return the total intramolecular energy
c     of residue ires in protein iprot, and will also return
c     the components in bndnrg, angnrg, dihnrg, ubynrg and nbnrg...
c
c     (bearing in mind whether the intra energy of any bonded residues
c     had already been calculated!)
c
c     (C) Christopher Woods
c
c########################################################################

      integer iprot,ires,inres
      double precision bndnrg,angnrg,dihnrg,ubynrg,nbnrg
      double precision bndnrgf,angnrgf,dihnrgf,ubynrgf,nbcnrgf,nbljnrgf
      double precision bndnrgb,angnrgb,dihnrgb,ubynrgb,nbcnrgb,nbljnrgb
      double precision ibndf,iangf,idihf,iubyf,ibndb,iangb,idihb,iubyb
      double precision icnrgf,icnrgb,iljnrgf,iljnrgb,imindist
      double precision nbcnrg,nbljnrg,icnrg,iljnrg
      integer i,inbnds,inangs,indihs,inubys,ibnd,iang,idih,iuby
      double precision bondEnergy,angleEnergy,dihedralEnergy,residueResidueEnergy
      double precision residueIntraNBEnergy,ureyBradleyEnergy

      residueIntraEnergy = ZERO
      bndnrg = ZERO
      angnrg = ZERO
      dihnrg = ZERO
      ubynrg = ZERO
      nbnrg = ZERO
      nbcnrg = ZERO
      nbljnrg = ZERO
      bndnrgf = ZERO
      bndnrgb = ZERO
      angnrgb = ZERO
      angnrgf = ZERO
      dihnrgb = ZERO
      dihnrgf = ZERO
      ubynrgb = ZERO
      ubynrgf = ZERO
      nbcnrgb = ZERO
      nbcnrgf = ZERO
      nbljnrgf = ZERO
      nbljnrgb = ZERO

      if (iprot.le.0 .or. iprot.gt.NProteins .or. ires.le.0) return
      if (ires.gt.NResidues(iprot)) return
      
c     if the energy of this residue has already been calculated, then return  
      if (DoneResNrg(iprot,ires)) return
      
c     run over all of the bonds in the residue...
      inbnds = NResBonds(iprot,ires)
      if (inbnds.ne.0) then
        do i=1,inbnds
c         make sure that we are not double adding the bond...
          ibnd = ResidueBond(iprot,ires,i)
          if (    .not.DoneResNrg(iprot,Bond(iprot,ibnd,4))
     .     .and.  .not.DoneResNrg(iprot,Bond(iprot,ibnd,5))) then
              bndnrg = bndnrg + bondEnergy(PROTEIN,iprot,ibnd,ibndb,ibndf)
              bndnrgf = bndnrgf + ibndf
              bndnrgb = bndnrgb + ibndb
          endif
        enddo
      endif
      
c     run over all of the angles in the residue...
      inangs = NResAngles(iprot,ires)
      if (inangs.ne.0) then
        do i=1,inangs
          iang = ResidueAngle(iprot,ires,i)
c         make sure that we don't double add angles...
          if (   .not.DoneResNrg(iprot,Angle(iprot,iang,5))
     .     .and. .not.DoneResNrg(iprot,Angle(iprot,iang,6))
     .     .and. .not.DoneResNrg(iprot,Angle(iprot,iang,7))) then
              angnrg = angnrg + angleEnergy(PROTEIN,iprot,iang,iangb,iangf)
              angnrgf = angnrgf + iangf
              angnrgb = angnrgb + iangb
          endif
        enddo
      endif
      
c     run over all of the dihedrals in the residue...
      indihs = NResDihedrals(iprot,ires)
      if (indihs.ne.0) then
        do i=1,indihs
          idih = ResidueDihedral(iprot,ires,i)
          if (   .not.DoneResNrg(iprot,Dihedral(iprot,idih,6))
     .     .and. .not.DoneResNrg(iprot,Dihedral(iprot,idih,7))
     .     .and. .not.DoneResNrg(iprot,Dihedral(iprot,idih,8))
     .     .and. .not.DoneResNrg(iprot,Dihedral(iprot,idih,9))) then
              dihnrg = dihnrg + dihedralEnergy(PROTEIN,iprot,idih,idihb,idihf)
              dihnrgb = dihnrgb + idihb
              dihnrgf = dihnrgf + idihf
          endif
        enddo
      endif
      
c     run over all of the Urey-Bradley terms in the residue...
      inubys = NResUreyBradleys(iprot,ires)
      if (inubys.ne.0) then
        do i=1,inubys
c         make sure that we are not double adding the uby...
          iuby = ResidueUreyBradley(iprot,ires,i)
          if (    .not.DoneResNrg(iprot,UreyBradley(iprot,iuby,4))
     .     .and.  .not.DoneResNrg(iprot,UreyBradley(iprot,iuby,5))) then
              ubynrg = ubynrg + ureyBradleyEnergy(PROTEIN,iprot,iuby,iubyb,iubyf)
              ubynrgf = ubynrgf + iubyf
              ubynrgb = ubynrgb + iubyb
          endif
        enddo
      endif

c     now get the intraresidue non-bonded term
      nbnrg = residueIntraNBEnergy(iprot,ires,nbcnrg,nbljnrg,nbcnrgb,nbljnrgb,
     .                                                     nbcnrgf,nbljnrgf)

c     now run over all residues, and calculate the res-res energy
      inres = NResidues(iprot)
!$OMP PARALLEL DO DEFAULT(NONE) ORDERED
!$OMP&SHARED(ires,iprot,inres)
!$OMP&PRIVATE(icnrg,iljnrg,icnrgf,icnrgb,iljnrgf,iljnrgb,imindist)
!$OMP&REDUCTION(+:nbnrg,nbcnrg,nbljnrg,nbcnrgf,nbcnrgb,nbljnrgf,nbljnrgb)
      do i=1,inres
         if (i.ne.ires) then
!$OMP ORDERED
           nbnrg = nbnrg + residueResidueEnergy(iprot,ires,i,
     .                     icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
!$OMP END ORDERED
           nbcnrg = nbcnrg + icnrg
           nbljnrg = nbljnrg + iljnrg
           nbcnrgf = nbcnrgf + icnrgf
           nbljnrgf = nbljnrgf + iljnrgf
           nbcnrgb = nbcnrgb + icnrgb
           nbljnrgb = nbljnrgb + iljnrgb
         endif
      enddo
!$OMP END PARALLEL DO

10    format('Energies= ',3f8.3)
      write(printstring,10) nbnrg,nbcnrg,nbljnrg
      call printLine(RESULTS,printstring)

      residueIntraEnergy = bndnrg + angnrg + dihnrg + ubynrg + nbnrg
      
      return 
      
      end
