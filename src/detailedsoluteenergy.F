      subroutine detailedSoluteEnergy(isol)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'enums.inc'
c###############################################################################
c
c       This subroutine calculates the same energy as 'soluteEnergy', but
c       instead prints out *all* of the energy components. This is useful
c       for template debugging.
c
c       (C) Christopher Woods  23 May 2003
c
c################################################################################

      integer i,isol
      integer savenrg
      integer lgth,lineLength
      double precision solenergy
      double precision soluteDetailedIntraEnergy,soluteSoluteEnergy
      double precision proteinSoluteEnergy,soluteSolventEnergy
      double precision icnrg,iljnrg
      double precision imindist
      double precision inrg
      double precision icnrgf,icnrgb,iljnrgf,iljnrgb

      savenrg = DONTCARE
      solenergy = ZERO

      if (isol.le.0 .or. isol.gt.NSolutes) then
        write(printstring,*) "Cannot print energy of invalid solute ",isol
        call printLine(INFO,printstring)
        call printLine(WARNING,printstring)
        return
      endif

      lgth = lineLength(SoluteName(isol))
      write(printstring,*) "Detailed energy of solute ",isol," ",SoluteName(isol)(1:lgth)
      call printLine(ENERGY,printstring)
            
c     the solute intramolecular energy
      inrg = soluteDetailedIntraEnergy(isol)
      write(printstring,*) "Solute intra-energy = ",inrg
      call printLine(ENERGY,printstring)
      solEnergy = solEnergy + inrg
      
c     inter-solute energy
      if (NSolutes.gt.1) then
        do i=1,NSolutes
          if (i.ne.isol .and. .not.IgnoreSolute(i)) then
            inrg =  soluteSoluteEnergy(isol,i,
     .             icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
            write(printstring,*) "Solute-solute ",i," energy = ",inrg
            call printLine(ENERGY,printstring)
            write(printstring,*) "Coulomb = ",icnrg," LJ = ",iljnrg
            call printLine(ENERGY,printstring)
            solEnergy = solEnergy + inrg
          endif
        enddo
      endif
      
c     solute-protein energy
      if (NProteins.gt.0) then
        do i=1,NProteins
          inrg = proteinSoluteEnergy(i,isol,
     .              icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
          write(printstring,*) "Solute-Protein ",i," energy = ",inrg
          call printLine(ENERGY,printstring)
          write(printstring,*) "Coulomb = ",icnrg," LJ = ",iljnrg
          call printLine(ENERGY,printstring)
          solEnergy = solEnergy + inrg
        enddo
      endif
      
c     solute-solvent energy
      if (NSolvents.gt.0) then
        inrg =  soluteSolventEnergy(savenrg,isol,
     .               icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist)
      
        write(printstring,*) "Solute-solvent energy = ",inrg
        call printLine(ENERGY,printstring)
        write(printstring,*) "Coulomb = ",icnrg," LJ = ",iljnrg
        call printLine(ENERGY,printstring) 
        write(printstring,*) "Mindist = ",sqrt(imindist)
        call printLine(ENERGY,printstring)
        solEnergy = solEnergy + inrg
      endif

      return
      end
