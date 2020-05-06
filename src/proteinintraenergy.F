      double precision function proteinIntraEnergy(iprot,
     .                      bndnrg, angnrg, dihnrg, ubynrg, nbcnrg, nbljnrg,
     .                      bndnrgb,angnrgb,dihnrgb,ubynrgb,nbcnrgb,nbljnrgb,
     .                      bndnrgf,angnrgf,dihnrgf,ubynrgf,nbcnrgf,nbljnrgf)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'constants.inc'
      include 'proteins.inc'
c#######################################################################
c
c     This function returns the total intramolecular energy of 
c     protein iprot. 
c
c     It also returns the components of the nrg, as well as the 
c     forwards and backwards FEP energies
c
c     *NOTE* That there are no checks that the energies of residues
c     within the protein have not been done before (since this is a 
c     'whole' molecule intramolecular calculation)
c
c     In other words, don't call this function if you are also calling
c     residueEnergy(iprot,ires) on iprot!
c
c     (C) Christopher Woods   7 February 2002
c
c#######################################################################

      integer iprot
      double precision bndnrg,angnrg,dihnrg,ubynrg,nbcnrg,nbljnrg
      double precision bndnrgb,angnrgb,dihnrgb,ubynrgb,nbcnrgb,nbljnrgb
      double precision bndnrgf,angnrgf,dihnrgf,ubynrgf,nbcnrgf,nbljnrgf
      double precision ibndb,ibndf,iangb,iangf,idihb,idihf,iubyb,iubyf
      double precision icnrg,iljnrg,icnrgf,icnrgb,iljnrgf,iljnrgb
      double precision residueResidueEnergy,imindist
      integer i,j,nres,nbnds,nangs,ndihs,nubys
      double precision bondEnergy,angleEnergy,dihedralEnergy,ureyBradleyEnergy
      double precision residueIntraNBEnergy

      proteinIntraEnergy = ZERO     

      bndnrg = ZERO
      angnrg = ZERO
      dihnrg = ZERO
      ubynrg = ZERO
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
      
      if (iprot.le.0 .or. iprot.gt.NProteins) return

      nres = NResidues(iprot)
      if (nres.le.0) return

c     the internal angles should have been updated when the protein
c     was built, or when it was moved!

c     loop over all of the bonds of this protein
      nbnds = NBonds(iprot)
      if (nbnds.gt.0) then
        do i=1,nbnds
          bndnrg = bndnrg + bondEnergy(PROTEIN,iprot,i,ibndb,ibndf)
          bndnrgb = bndnrgb + ibndb
          bndnrgf = bndnrgf + ibndf
        enddo
      endif
      
      nangs = NAngles(iprot)
      if (nangs.gt.0) then
        do i=1,nangs
          angnrg = angnrg + angleEnergy(PROTEIN,iprot,i,iangb,iangf)
          angnrgb = angnrgb + iangb
          angnrgf = angnrgf + iangf
        enddo
      endif
     
      ndihs = NDihedrals(iprot)
      if (ndihs.gt.0) then
        do i=1,ndihs
          dihnrg = dihnrg + dihedralEnergy(PROTEIN,iprot,i,idihb,idihf)
          dihnrgb = dihnrgb + idihb
          dihnrgf = dihnrgf + idihf
        enddo
      endif
        
c     loop over all of the urey-bradley terms of this protein
      nubys = NUreyBradleys(iprot)
      if (nubys.gt.0) then
        do i=1,nubys
          ubynrg = ubynrg + ureyBradleyEnergy(PROTEIN,iprot,i,iubyb,iubyf)
          ubynrgb = ubynrgb + iubyb
          ubynrgf = ubynrgf + iubyf
        enddo
      endif

c     now run over all residue pairs in the protein and calculate the
c     non-bonded energy... - this allows residuebased cutting
      nres = NResidues(iprot)
      if (nres.gt.1) then
        do i=1,nres-1
          do j=i+1,nres
c           don't bother if this is between a pair of fixed residues
            if (.not.(FixedRes(iprot,i).and.FixedRes(iprot,j))) then
              proteinIntraEnergy = proteinIntraEnergy + residueResidueEnergy(iprot,i,j,
     .                  icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
              nbcnrg = nbcnrg + icnrg
              nbljnrg = nbljnrg + iljnrg
              nbcnrgf = nbcnrgf + icnrgf
              nbcnrgb = nbcnrgb + icnrgb
              nbljnrgf = nbljnrgf + iljnrgf
              nbljnrgb = nbljnrgb + iljnrgb
            endif
          enddo
        enddo
      endif
      
c     this, however, has missed out the intra-residue nonbonded interactions
c     loop over each residue and add those on
      do i=1,nres
c       don't bother if this is a fixed residue
        if (.not.FixedRes(iprot,i)) then
          proteinIntraEnergy = proteinIntraEnergy + residueIntraNBEnergy(iprot,i,
     .                   icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf)
          nbcnrg = nbcnrg + icnrg
          nbljnrg = nbljnrg + iljnrg
          nbcnrgf = nbcnrgf + icnrgf
          nbcnrgb = nbcnrgb + icnrgb
          nbljnrgf = nbljnrgf + iljnrgf
          nbljnrgb = nbljnrgb + iljnrgb
        endif
      enddo

      proteinIntraEnergy = bndnrg + angnrg + dihnrg + ubynrg + nbcnrg + nbljnrg
      
      return
      end
