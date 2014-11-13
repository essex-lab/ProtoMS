      
      double precision function extraEnergy(savenrg,nrgb,nrgf)
      include 'dimensions.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'constants.inc'
      include 'printer.inc'
      include 'extraenergies.inc'
      include 'solutes.inc'
c###########################################################
c
c     This is the master function responsible for calculating
c     all of the extra energy terms (restraints, umbrellas,
c     etc.). It returns the total 'extraEnergy' for the 
c     current system configuration, and the backwards
c     and forwards energies. If savenrg=OLD then the 
c     energy is saved as an 'old' energy, if savenrg=NEW
c     then saved as 'new' energy, else lost to the ether!
c
c     Note that this function will only calculate the extra
c     energies of parts of the system that have moved.
c     You will need to call 'setAllMoved' if you want the
c     total extra energy.
c
c     (C) Christopher Woods, November 2004
c
c###########################################################

      integer savenrg
      double precision TotalSurface
      double precision temp,surfacef,surfaceb
      double precision sanrg,sanrgf,sanrgb
      double precision GBEnergy,gbnrgf,gbnrgb
      double precision nrgb,nrgf,inrgf,inrgb
      double precision solventCapEnergy,harmonicRestraintEnergy
      double precision flatbottomRestraintEnergy
      double precision DihedralHarmonicRestraintEnergy
      double precision nrg
      integer i
      logical dualsasa
      double precision BndharmonicRestraintEnergy,hardwallEnergy

      extraEnergy = ZERO
      nrgb = ZERO
      nrgf = ZERO

      inrgf = ZERO
      inrgb = ZERO
      temp = ZERO
    
c     go through each type of extra energy term and see whether or not it is activated...
      if (SolventCap) then
c      calculate the solvent cap energy
       extraEnergy = extraEnergy + solventCapEnergy(savenrg,inrgb,inrgf)
       nrgb = nrgb + inrgb
       nrgf = nrgf + inrgf
c       write(printstring,*) 'The svn cap Energy ',nrgb
c       call printLine(DEBUG,printstring)
      endif            

      if (NHarmonicRestraints.gt.0) then
c       calculate the energy of harmonic restraints
        extraEnergy = extraEnergy + harmonicRestraintEnergy(savenrg,inrgb,inrgf)
        nrgb = nrgb + inrgb
        nrgf = nrgf + inrgf
      endif

      if (NFlatBottomRestraints.gt.0) then
c       calculate the energy of flat bottom harmonic restraints
        extraEnergy = extraEnergy + flatbottomRestraintEnergy(savenrg,inrgb,inrgf)
        nrgb = nrgb + inrgb
        nrgf = nrgf + inrgf
c        write(printstring,*) 'FBOT NRG ',nrgb
c        call printLine(DEBUG,printstring)
      endif

      if (NDihedralHarmonicRestraints.gt.0) then
c       extraEnergy = extraEnergy + DihedralHarmonicRestraintEnergy(savenrg,inrgb,inrgf)
        nrg = DihedralHarmonicRestraintEnergy(savenrg,inrgb,inrgf)
        extraEnergy = extraEnergy + nrg

        nrgb = nrgb + inrgb
        nrgf = nrgf + inrgf
      endif

      if (nBndHarmonicRestraints .gt. 0) then 
        nrg = BndHarmonicRestraintEnergy(savenrg,inrgb,inrgf)
        extraEnergy = extraEnergy + nrg
        nrgb = nrgb + inrgb
        nrgf = nrgf + inrgf        
      endif

      if (NHardwallPotential.gt.0) then
c       calculate the energy of hardwall restraints
        extraEnergy = extraEnergy + hardwallEnergy(savenrg,inrgb,inrgf)
        nrgb = nrgb + inrgb
        nrgf = nrgf + inrgf
      endif


c     if (NUmbrellas.gt.0) then
c       extraEnergy = extraEnergy + umbrellaEnergy(savenrg,inrgb,inrgf)
c       nrgb = nrgb + inrgb
c       nrgf = nrgf + inrgf
c     endif

      if (doSA) then
c       Calculate the Surface Area of the System
c       JM june 06. Coding SASA support for dual topology
c       First detect if there is a pair of dual topology solutes
        dualsasa = .false.
        if (NSolutes .gt. ZERO) then
          do 200 i=1,NSolutes
            if (IgnoreSolute(i)) goto 200
            if (DualTopologySolute(i).gt.0) then
              call dualtopologysasa(i,savenrg,sanrg,sanrgf,sanrgb,temp,surfacef,surfaceb)
              dualsasa = .true.
c             Break out the loop. THIS WILL BUG SIMULATIONS WITH MORE THAN ONE PAIR OF DUAL TOPOLOGY SOLUTES 
c              STOP
              goto 210
            endif
 200      continue
 210      continue
        endif        
        if (.not. dualsasa) temp = TotalSurface(savenrg,surfacef,surfaceb,sanrg,sanrgf,sanrgb)
c        if (savenrg .eq. OLD) then 
c          write(printstring,*) ' OLD SASA adding ',sanrg,' ',sanrgf,' ',sanrgb,' surf ',temp,' F ',surfacef,' B ',surfaceb
c        else if (savenrg .eq. NEW) then 
c          write(printstring,*) ' NEW SASA adding ',sanrg,' ',sanrgf,' ',sanrgb,' surf ',temp,' F ',surfacef,' B ',surfaceb
c        endif
c        call printLine(DEBUG,printstring)           
        extraEnergy = extraEnergy + sanrg
        nrgb = nrgb + sanrgb
        nrgf = nrgf + sanrgf
      endif  

      if (doGB) then
        temp = GBEnergy(savenrg,gbnrgf,gbnrgb)
c        if (savenrg .eq. OLD) then 
c          write(printstring,*) ' ## Old GB nrg ',temp,' ',gbnrgf,' ',gbnrgb
c        else if (savenrg .eq. NEW) then
c          write(printstring,*) ' ## New GB nrg ',temp,' ',gbnrgf,' ',gbnrgb
c        endif
c        call printLine(DEBUG,printstring)
        extraEnergy = extraEnergy + temp
        nrgb = nrgb + gbnrgb
        nrgf = nrgf + gbnrgf
      endif

c     see whether or not we need to save the energies
      if (savenrg.eq.OLD) then
        oldExtraEnergyPart = extraEnergy
c        write(printstring,*) 'OLD extraNrg ',extraEnergy
c        call printLine(DEBUG,printstring)
        oldExtraEnergyPartB = nrgb
        oldExtraEnergyPartF = nrgf
      else if (savenrg.eq.NEW) then
        newExtraEnergyPart = extraEnergy
c        write(printstring,*) 'NEW extraNrg ',extraEnergy
c        call printLine(DEBUG,printstring)
        newExtraEnergyPartB = nrgb
        newExtraEnergyPartF = nrgf
      endif
      
      return
      
      end
