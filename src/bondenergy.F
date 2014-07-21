      double precision function bondEnergy(ityp,ival,ibnd,bondEnergyB,bondEnergyF)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'solutes.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'parameters.inc'
c#######################################################
c
c     This is the function used to evaluate the energy
c     of bond ibnd in ityp ic. It returns the 
c     energy in kcal mol-1, or ZERO if the bond
c     doesn't exist...
c
c     (C) Christopher Woods   23 January 2002
c
c#######################################################

      integer iprot,ibnd,isol,ityp,iparam,ival,ires1,ires2
      double precision ir,ik,bondEnergyB,bondEnergyF,ibndval
      double precision ibndvalf,ibndvalb,irf,ikf,irb,ikb
      logical ivariable
      
c     BndEnergy uses harmonic potential...
c     BndEnergy = 0.5 * k * (R - R0)**2      
                  
      BondEnergy = ZERO
      BondEnergyF = ZERO
      BondEnergyB = ZERO                  
                       
      ivariable = .false.
      if (ityp.eq.PROTEIN) then
        iprot = ival
        ires1 = Zid2ResID(iprot,Bond(iprot,ibnd,1),1)
        ires2 = Zid2ResID(iprot,Bond(iprot,ibnd,2),1)
c       don't bother evaluating bonds which are within fixed residues   
        if ( FixedRes(iprot,ires1) .and. FixedRes(iprot,ires2) ) then
           return
        endif
        iparam = Bond(iprot,ibnd,3)
        ibndval = BondLength(iprot,ibnd)
        ibndvalf = ibndval
        ibndvalb = ibndval
      else if (ityp.eq.SOLUTE) then
        isol = ival
        iparam = SolBond(isol,ibnd,3)
        ibndval = SolBondLength(isol,ibnd)
        if (NSolVariable(isol).gt.0.and.DifferentLambdas) then
          ivariable = .true.
          ibndvalf = SolBondLength(NSolForwards(isol),ibnd)
          ibndvalb = SolBondLength(NSolBackwards(isol),ibnd)
        else
          ibndvalf = ibndval
          ibndvalb = ibndval
        endif
      else 
        isol = 0
        iprot = 0
        iparam = 0
        ibndval = ZERO
        ibndvalf = ibndval
        ibndvalb = ibndval
        write(printstring,*) "Unknown bondEnergy() system type requested ",ityp
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is probably a programming bug!")
        return
      endif
      
      if (iparam.lt.0) then
        ik = VariableBndParam(-iparam,2,1)
        ir = VariableBndParam(-iparam,2,2)
      else if (iparam.eq.0) then
        bondEnergy = ZERO
        bondEnergyF = ZERO
        bondEnergyB = ZERO
        return
      else 
        ik = BndParam(iparam,1)
        ir = BndParam(iparam,2)
      endif

      bondEnergy = ik * (ibndval-ir)**2

      if (DifferentLambdas) then
        if (iparam.gt.0 .and. .not.ivariable) then
          bondEnergyf = bondEnergy
          bondEnergyb = bondEnergy
        else
          if (iparam.lt.0) then
c           backwards
            ikb = VariableBndParam(-iparam,4,1)
            irb = VariableBndParam(-iparam,4,2)
c           forwards
            ikf = VariableBndParam(-iparam,5,1)
            irf = VariableBndParam(-iparam,5,2)
          else
            ikf = ik
            ikb = ik
            irf = ir
            irb = ir
          endif
          bondEnergyB = ikb * (ibndvalb-irb)**2
          bondEnergyF = ikf * (ibndvalf-irf)**2
        endif
      else
        bondenergyF = bondEnergy
        bondEnergyB = bondEnergy
      endif 
      
      return
      
      end
