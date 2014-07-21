      subroutine detailedResidueEnergy(iprot,ires)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'simparams.inc'
      include 'printer.inc'
c#############################################################
c
c     This is identical to residueEnergy, except that
c     this prints out all of the energy components
c     (much like detailedSoluteEnergy)
c
c     (C) Christopher Woods, November 2004
c
c#############################################################

      integer iprot,ires,i
      integer SAVENRG
      double precision inrg,resnrg,mindist
      double precision icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf
      double precision getLowScale,getHighScale,scl
      double precision residueDetailedIntraEnergy,residueProteinEnergy
      double precision residueSoluteEnergy,residueSolventEnergy
      
      SAVENRG = DONTCARE
      resnrg = ZERO
      inrg = ZERO
      
      if (iprot.le.0 .or. iprot.gt.NProteins .or. ires.le.0) return
      if (ires.gt.NResidues(iprot)) return
      
      inrg = residueDetailedIntraEnergy(iprot,ires)
      
      write(printstring,*) "Residue intra energy = ",inrg
      call prettyLine(ENERGY,printstring)
      
      resnrg = resnrg + inrg
      
c     residue-protein energies
      if (NProteins.gt.1) then
        do i=1,NProteins
          if (i.ne.iprot) then
            inrg = residueProteinEnergy(savenrg,iprot,ires,i,icnrg,iljnrg,icnrgb,iljnrgb,
     .                                  icnrgf,iljnrgf,mindist)
            write(printstring,*) "Residue - protein ",i," energy = ",inrg
            call printLine(ENERGY,printstring)
            write(printstring,*) "Coulomb = ",icnrg," LJ = ",iljnrg
            call printLine(ENERGY,printstring)
            resnrg = resnrg + inrg
          endif
        enddo
      endif
      
c     residue-solute energies
      if (NSolutes.gt.0) then
        do i=1,NSolutes
          inrg = residueSoluteEnergy(iprot,ires,i,
     .                icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,mindist)
c         see if this is a dual topology solute, and thus whether we need to
c         scale the energies
          if (DualTopologySolute(i).gt.0) then
c           this is the low solute
            scl = getLowScale(Lambda)
c            sclb = getLowScale(LambdaB)
c            sclf = getLowScale(LambdaF)
          else if (DualTopologySolute(i).eq.0) then
c           this is not a dual topology solute
            scl = ONE
c            sclb = ONE
c            sclf = ONE
          else
c           this is the high solute
            scl = getHighScale(Lambda)
c            sclb = getHighScale(LambdaB)
c            sclf = getHighScale(LambdaF)
          endif           

          inrg = inrg * scl
          write(printstring,*) "Residue - solute ",i," energy = ",inrg
          call printLine(ENERGY,printstring)
          write(printstring,*) "Coulomb = ",icnrg*scl," LJ = ",iljnrg*scl
          call printLine(ENERGY,printstring)
          resnrg = resnrg + inrg
        enddo
      endif
      
c     residue-solvent energy
      if (NSolvents.gt.0) then
        inrg = residueSolventEnergy(savenrg,iprot,ires,
     .           icnrg,iljnrg,icnrgb,iljnrgb,icnrgf,iljnrgf,mindist)
        
        write(printstring,*) "Residue - solvent energy = ",inrg
        call printLine(ENERGY,printstring)
        write(printstring,*) "Coulomb = ",icnrg," LJ = ",iljnrg
        call printLine(ENERGY,printstring)
      
        resnrg = resnrg + inrg
      endif
      
      write(printstring,*) "Total residue energy = ",resnrg
      call prettyLine(ENERGY,printstring)
      
      return
      end
