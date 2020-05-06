      subroutine detailedSoluteEnergy(isol)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'solutes.inc'
      include 'simparams.inc'
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
c       Modified by Samuel Genheden, 2012-2013
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
      double precision scl,sclf,sclb,getHighScale,getLowScale

      savenrg = DONTCARE
      solenergy = ZERO


      if (isol.le.0 .or. isol.gt.NSolutes) then
        write(printstring,*) "Cannot print energy of invalid solute ",isol
        call printLine(INFO,printstring)
        call printLine(WARNING,printstring)
        return
      endif

c     see if this is a dual topology solute - get the scale factor if it is
      if (DualTopologySolute(isol).eq.0) then
c       not a dual topology solute
        scl = ONE
        sclf = ONE
        sclb = ONE
      else if (DualTopologySolute(isol).gt.0) then
c       this is the dual topology solute at lambda=0
        scl = getLowScale(Lambda)
        sclb = getLowScale(LambdaB)
        sclf = getLowScale(LambdaF)
      else
c       this is the dual topology solute at lambda=1
        scl = getHighScale(Lambda)
        sclb = getHighScale(LambdaB)
        sclf = getHighScale(LambdaF)
      endif  
      scl = scl*SoluteTheta(isol) ! SG
      sclf = sclf*SoluteTheta(isol) ! SG
      sclb = sclb*SoluteTheta(isol) ! SG

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
!       SG Added DONTCARE
          inrg = proteinSoluteEnergy(savenrg,i,isol,
     .              icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist,scl,sclf,sclb)
c         If the soft-core is on, the coulombic is scaled further here
          if (Solutesoftcore(isol)) then
            icnrg = icnrg * (scl)**CoulSoftPower
          endif                
          icnrg = icnrg*scl
          iljnrg =  iljnrg*scl
!         SG
          inrg = icnrg + iljnrg
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
     .               icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,imindist,scl,sclf,sclb)
c       If the soft-core is on, the coulombic is scaled further here
        if (Solutesoftcore(isol)) then
          icnrg = icnrg * (scl)**CoulSoftPower
        endif                
        icnrg = icnrg*scl
        iljnrg =  iljnrg*scl
!       SG
        inrg = icnrg + iljnrg
        write(printstring,*) "Solute-solvent energy = ",inrg
        call printLine(ENERGY,printstring)
        write(printstring,*) "Coulomb = ",icnrg," LJ = ",iljnrg
        call printLine(ENERGY,printstring)
        solEnergy = solEnergy + inrg
      endif

      return
      end
