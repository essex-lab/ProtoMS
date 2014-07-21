      double precision function soluteDetailedIntraEnergy(isol)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'simparams.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'enums.inc'
      include 'parameters.inc'
      include 'flags.inc'
c#######################################################################
c
c     This function returns the total intramolecular energy of 
c     solute isol, printing the detailed list of components
c
c     (C) Christopher Woods   23 May 2003
c
c#######################################################################

      integer isol
      double precision bndnrg,angnrg,dihnrg,ubynrg,nbcnrg,nbljnrg,mindist
      double precision bndnrgf,angnrgf,dihnrgf,ubynrgf,nbcnrgf,nbljnrgf
      double precision bndnrgb,angnrgb,dihnrgb,ubynrgb,nbcnrgb,nbljnrgb
      double precision ibndb,ibndf,iangb,iangf,idihb,idihf,iubyb,iubyf
      double precision bondEnergy,angleEnergy,dihedralEnergy,ureyBradleyEnergy
      double precision inrg
      double precision rad2Deg
      integer inbnds,inangs,indihs,inubys,i,innb,ic
      include 'nbenergyvar2.inc'
      double precision charge,sigma,epslon      
      character*10 inam1,inam2,inam3,inam4,getSolAtomName
c     temporary cheat
      double precision scl,sclf,sclb
      scl = 1.0d0
      sclf = 1.0d0
      sclb = 1.0d0
c     end of cheat
      soluteDetailedIntraEnergy = ZERO     
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
c     calculate the bond / angle / dihedral energies
      inbnds = NSolBonds(isol)
      if (inbnds.gt.0) then
        ic = 0
        do i=1,inbnds
          inrg = bondEnergy(SOLUTE,isol,i,ibndb,ibndf)
          if (SolBond(isol,i,3).ne.0) then
            ic = ic + 1
            inam1 = getSolAtomName(isol,SolBond(isol,i,1))
            inam2 = getSolAtomName(isol,SolBond(isol,i,2))
            
            write(printstring,10) inam1,inam2,
     .                            SolBondLength(isol,i),inrg
10          format('Bond ',a,' - ',a,' Length = ',f10.5,' Energy = ',f10.5)
            call printLine(ENERGY,printstring)
          endif
          bndnrg = bndnrg + inrg
        enddo
      
c       print out the total energy if there were any energetic bonds        
        if (ic.gt.0) then
          write(printstring,*) "Total bond energy = ",bndnrg
          call prettyLine(ENERGY,printstring)
        endif
      endif
      
      inangs = NSolAngles(isol)
      if (inangs.gt.0) then
        ic = 0
        do i=1,inangs
          inrg = angleEnergy(SOLUTE,isol,i,iangb,iangf)
          if (SolAngle(isol,i,4).ne.0) then
            ic = ic + 1
            inam1 = getSolAtomName(isol,SolAngle(isol,i,1))
            inam2 = getSolAtomName(isol,SolAngle(isol,i,2))
            inam3 = getSolAtomName(isol,SolAngle(isol,i,3))
            
            write(printstring,11) inam1,inam2,inam3,
     .                            rad2Deg(SolAngleSize(isol,i)),inrg
11          format('Angle ',a,' - ',a,' - ',a,' Size = ',f10.5,' Energy = ',f10.5)
            call printLine(ENERGY,printstring)
          endif
          angnrg = angnrg + inrg
        enddo
      
        if (ic.gt.0) then
          write(printstring,*) "Total angle energy = ",angnrg
          call prettyLine(ENERGY,printstring)
        endif
      endif
      
      indihs = NSolDihedrals(isol)
      if (indihs.gt.0) then
        ic = 0
        do i=1,indihs
          inrg = dihedralEnergy(SOLUTE,isol,i,idihb,idihf)
          if (SolDihedral(isol,i,5).ne.0) then
            ic = ic + 1
            inam1 = getSolAtomName(isol,SolDihedral(isol,i,1))
            inam2 = getSolAtomName(isol,SolDihedral(isol,i,2))
            inam3 = getSolAtomName(isol,SolDihedral(isol,i,3))
            inam4 = getSolAtomName(isol,SolDihedral(isol,i,4))
            
            write(printstring,12) inam1,inam2,inam3,inam4,
     .                            rad2Deg(SolDihedralSize(isol,i)),inrg
12          format('Dihedral ',a ,' - ',a,' - ',a,' - ',a,' Size = ',f10.5,' Energy = ',f10.5)
            call printLine(ENERGY,printstring)
          endif
          dihnrg = dihnrg + inrg
        enddo
        
        if (ic.gt.0) then
          write(printstring,*) "Total dihedral energy = ",dihnrg
          call prettyLine(ENERGY,printstring)
        endif
      endif
      
      inubys = NSolUreyBradleys(isol)
      if (inubys.gt.0) then
        ic = 0
        do i=1,inubys
          inrg = ureyBradleyEnergy(SOLUTE,isol,i,iubyb,iubyf)
          if (SolUreyBradley(isol,i,4).ne.0) then
            ic = ic + 1
            inam1 = getSolAtomName(isol,SolUreyBradley(isol,i,1))
            inam2 = getSolAtomName(isol,SolUreyBradley(isol,i,2))
            inam3 = getSolAtomName(isol,SolUreyBradley(isol,i,3))
            
            write(printstring,13) inam1,inam2,inam3,
     .                            SolUreyBradleySize(isol,i),inrg
13          format('Urey-Bradley ',a,' - ',a,' - ',a,' Size = ',f10.5,' Energy = ',f10.5)
            call printLine(ENERGY,printstring)
          endif
          ubynrg = ubynrg + inrg
        enddo
        
        if (ic.gt.0) then
          write(printstring,*) "Total Urey-Bradley energy = ",ubynrg
          call prettyLine(ENERGY,printstring)
        endif
      endif
      
c     now run over all of the non-bonded pairs in the solute molecule
      innb = NSolNB(isol)

c     no periodic boundaries in an intramolecular NB calculation
      boundx = .false.
      boundy = .false.
      boundz = .false.
            
      if (innb.gt.0) then
        do i=1,innb
          funcisol1 = isol
          funcisol2 = isol
          funcat1 = SolNB(isol,i,1)
          funcat2 = SolNB(isol,i,2)
          inam1 = getSolAtomName(isol,funcat1)
          inam2 = getSolAtomName(isol,funcat2)
          
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
            write(printstring,14) 'NB14',inam1,inam2,
     .            dsqrt(funcmindist),OneOverScl14C*funccnrg,OneOverScl14LJ*funcljnrg,
     .              ' SCLC = ',OneOverScl14C,' SCLJ = ',OneOverScl14LJ
          else
            nbcnrg = nbcnrg + funccnrg
            nbljnrg = nbljnrg + funcljnrg
            nbcnrgf = nbcnrgf + funccnrgf
            nbcnrgb = nbcnrgb + funccnrgb
            nbljnrgf = nbljnrgf + funcljnrgf
            nbljnrgb = nbljnrgb + funcljnrgb
            write(printstring,14) 'NB  ',inam1,inam2,
     .                         dsqrt(funcmindist),funccnrg,funcljnrg
          endif
          
14        format(a4,1x,a," - ",a,' Distance = ',f10.5,' Coul = ',g15.5,' LJ = ',g15.5,
     .           a8,f8.3,a8,f8.3 )
c         if the distance is less than sqrt(6.0) A, then put a box around it
          if (funcmindist.lt.6.0) then
            call prettyLine(ENERGY,printstring)
          else          
            call printLine(ENERGY,printstring)
          endif
          
        enddo
        write(printstring,*) "Total nb energy = ",nbcnrg+nbljnrg," coul = (",nbcnrg,") lj = (",nbljnrg,")"
        call prettyLine(ENERGY,printstring)      
      endif
      
      soluteDetailedIntraEnergy = bndnrg + angnrg + dihnrg + ubynrg + nbcnrg + nbljnrg
      
      write(printstring,*) "Total solute intra energy = ",soluteDetailedIntraEnergy
      call prettyLine(ENERGY,printstring)
      
      return
      end
