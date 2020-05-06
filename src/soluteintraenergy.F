      double precision function soluteIntraEnergy(isol,
     .          bndnrg, angnrg, dihnrg, ubynrg, nbcnrg, nbljnrg,
     .          bndnrgb,angnrgb,dihnrgb,ubynrgb,nbcnrgb,nbljnrgb,
     .          bndnrgf,angnrgf,dihnrgf,ubynrgf,nbcnrgf,nbljnrgf)
      include 'dimensions.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'simparams.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'flags.inc'
      include 'quickpot.inc'
c#######################################################################
c
c     This function returns the total intramolecular energy of 
c     solute isol. 
c
c     (C) Christopher Woods   7 February 2002
c
c#######################################################################

      integer isol
      double precision bndnrg,angnrg,dihnrg,ubynrg,nbcnrg,nbljnrg,mindist
      double precision bndnrgf,angnrgf,dihnrgf,ubynrgf,nbcnrgf,nbljnrgf
      double precision bndnrgb,angnrgb,dihnrgb,ubynrgb,nbcnrgb,nbljnrgb
      double precision ibndb,ibndf,iangb,iangf,idihb,idihf,iubyb,iubyf
      double precision bondEnergy,angleEnergy,dihedralEnergy,ureyBradleyEnergy
      integer inbnds,inangs,indihs,inubys,i,innb
      include 'nbenergyvar2.inc'
      double precision charge,sigma,epslon      
c     The scale values are 1 since no softcore is applied to the intramolecular
c     energy
      double precision scl,sclf,sclb
      scl = 1.0d0
      sclf = 1.0d0
      sclb = 1.0d0
c     end of cheat
      soluteIntraEnergy = ZERO     
      bndnrg = ZERO
      angnrg = ZERO
      dihnrg = ZERO
      ubynrg = ZERO
      nbcnrg = ZERO
      nbljnrg = ZERO
      bndnrgf = ZERO
      angnrgf = ZERO
      dihnrgf = ZERO
      ubynrgf = ZERO
      nbcnrgf = ZERO
      nbljnrgf = ZERO
      bndnrgb = ZERO
      angnrgb = ZERO
      dihnrgb = ZERO
      ubynrgb = ZERO
      nbcnrgb = ZERO
      nbljnrgb = ZERO
      mindist = 99999.0
      softcore = .false.
      ddd = .false.
c     return if this solute has been done before
      if (DoneSolNrg(isol)) return

      if (IgnoreSolute(isol)) return
      
c     calculate the bond / angle / dihedral energies
      inbnds = NSolBonds(isol)
      if (inbnds.gt.0) then
        do i=1,inbnds
          bndnrg = bndnrg + bondEnergy(SOLUTE,isol,i,ibndb,ibndf)
          bndnrgf = bndnrgf + ibndf
          bndnrgb = bndnrgb + ibndb
        enddo
      endif
      if (doRest) then
        bndnrg = bndnrg*restscaling(1)
        bndnrgf = bndnrgf*restscaling(1)
        bndnrgb = bndnrgb*restscaling(1)
      endif
      
      inangs = NSolAngles(isol)
      if (inangs.gt.0) then
        do i=1,inangs
          angnrg = angnrg + angleEnergy(SOLUTE,isol,i,iangb,iangf)
          angnrgf = angnrgf + iangf
          angnrgb = angnrgb + iangb
        enddo
      endif
      if (doRest) then
        angnrg = angnrg*restscaling(1)
        angnrgf = angnrgf*restscaling(1)
        angnrgb = angnrgb*restscaling(1)
      endif
      
      indihs = NSolDihedrals(isol)
      if (indihs.gt.0) then
        do i=1,indihs
          dihnrg = dihnrg + dihedralEnergy(SOLUTE,isol,i,idihb,idihf)
          dihnrgf = dihnrgf + idihf
          dihnrgb = dihnrgb + idihb
        enddo
      endif
      if (doRest) then
        dihnrg = dihnrg*restscaling(2)
        dihnrgf = dihnrgf*restscaling(2)
        dihnrgb = dihnrgb*restscaling(2)
      endif

      inubys = NSolUreyBradleys(isol)
      if (inubys.gt.0) then
        do i=1,inubys
          ubynrg = ubynrg + ureyBradleyEnergy(SOLUTE,isol,i,iubyb,iubyf)
          ubynrgf = ubynrgf + iubyf
          ubynrgb = ubynrgb + iubyb
        enddo
      endif
      if (doRest) then
        ubynrg = ubynrg*restscaling(1)
        ubynrgf = ubynrgf*restscaling(1)
        ubynrgb = ubynrgb*restscaling(1)
      endif

c     now run over all of the non-bonded pairs in the solute molecule
      innb = NSolNB(isol)

c     no periodic boundaries in an intramolecular NB calculation
      boundx = .false.
      boundy = .false.
      boundz = .false.

c     Activate Distance Dependent Dielectric calculation
      if (QuickMove .and. DDDon) ddd = .true.

      if (innb.gt.0) then
        do i=1,innb
          funcisol1 = isol
          funcisol2 = isol
          funcat1 = SolNB(isol,i,1)
          funcat2 = SolNB(isol,i,2)
          
c         intialise the nb energy calculation
          include 'nbinitial.inc'
c         load up the solute-solute parameters
          include 'nbsolsol.inc'
c         actually perform the calculation
          pertchg = perturbeIntra
          include 'nbenergy2.inc'

          if (SolNB14(isol,i)) then
            nbcnrg = nbcnrg + OneOverScl14C*funccnrg
            nbljnrg = nbljnrg + OneOverScl14LJ*funcljnrg
            nbcnrgf = nbcnrgf + OneOverScl14C*funccnrgf
            nbcnrgb = nbcnrgb + OneOverScl14C*funccnrgb
            nbljnrgf = nbljnrgf + OneOverScl14LJ*funcljnrgf
            nbljnrgb = nbljnrgb + OneOverScl14LJ*funcljnrgb
          else
            nbcnrg = nbcnrg + funccnrg
            nbljnrg = nbljnrg + funcljnrg
            nbcnrgf = nbcnrgf + funccnrgf
            nbcnrgb = nbcnrgb + funccnrgb
            nbljnrgf = nbljnrgf + funcljnrgf
            nbljnrgb = nbljnrgb + funcljnrgb
          endif          
        enddo
      endif
      if (doREST) then
        nbcnrg = nbcnrg * restscaling(4)
        nbcnrgf = nbcnrgf * restscaling(4)
        nbcnrgb = nbcnrgb * restscaling(4)
        nbljnrg = nbljnrg * restscaling(3)
        nbljnrgf = nbljnrgf * restscaling(3)
        nbljnrgb = nbljnrgb * restscaling(3)
      endif
            
      soluteIntraEnergy = bndnrg + angnrg + dihnrg + ubynrg + nbcnrg + nbljnrg
c      write(printstring,*) 'bnd ',bndnrg,' ang ',angnrg,' dih ',dihnrg,
c     .  ' nbc ',nbcnrg,' nbl ',nbljnrg
c      call printLine(DEBUG,printstring)
c      write(printstring,*) 'Solute intra energy ',soluteIntraEnergy
c      call printLine(DEBUG,printstring)
c$$$      write(printstring,*) ' bndnrg ',bndnrg,' bndnrgF ',bndnrgF,
c$$$     &  ' bndnrgB ',bndnrgB,' angnrg ',angnrg,' angnrgF ',angnrgF,
c$$$     &  ' angnrgB ',angnrgB
c$$$      call printLine(DEBUG,printstring)
c$$$      write(printstring,*) ' dihnrg ',dihnrg,' dihnrgF ',dihnrgF,
c$$$     &  ' dihnrgB ',dihnrgB
c$$$      call printLine(DEBUG,printstring)
c$$$      write(printstring,*) ' nbljnrg ',nbljnrg,' nbljnrgF ',nbljnrgF,
c$$$     &  ' nbljnrgB ',nbljnrgB,' nbcnrg ',nbcnrg,' nbcnrgF ',nbcnrgF,
c$$$     &  ' nbcnrgB ',nbcnrgB
c$$$      call printLine(DEBUG,printstring)
      return
      end
