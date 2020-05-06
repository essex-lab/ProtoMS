      double precision function dihedralEnergy(ityp,ival,idih,dihedralEnergyB,dihedralEnergyF)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'solutes.inc'
      include 'simparams.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'parameters.inc'
c#######################################################
c
c     This is the function used to evaluate the energy
c     of dihedral idih in ityp ival. It returns the 
c     energy in kcal mol-1, or ZERO if the dihedral
c     doesn't exist...
c
c     (C) Christopher Woods   23 January 2002
c
c#######################################################

      integer iprot,idih,iparam,ityp,ival,isol,ires1,ires2,ires3,ires4
      integer iparb,iparf
      double precision ang
      double precision dihedralEnergyF,dihedralEnergyB
      double precision dihEnergy
      double precision angf,angb,nrg0,nrg1
      double precision perturbForceField

      DihedralEnergy = ZERO
      DihedralEnergyF = ZERO
      DihedralEnergyB = ZERO
                       
      if (ityp.eq.PROTEIN) then
        iprot = ival
c       don't bother if the dihedral is within fixed residues
        ires1 = Zid2ResID(iprot,Dihedral(iprot,idih,1),1)
        ires2 = Zid2ResID(iprot,Dihedral(iprot,idih,2),1)
        ires3 = Zid2ResID(iprot,Dihedral(iprot,idih,3),1)
        ires4 = Zid2ResID(iprot,Dihedral(iprot,idih,4),1)
        if (FixedRes(iprot,ires1).and.FixedRes(iprot,ires2).and.
     .      FixedRes(iprot,ires3).and.FixedRes(iprot,ires4)) then
          return
        endif
        iparam = Dihedral(iprot,idih,5)
        ang = DihedralSize(iprot,idih)
        angf = ang
        angb = ang
      else if (ityp.eq.SOLUTE) then
        isol = ival
        iparam = SolDihedral(isol,idih,5)
        ang = SolDihedralSize(isol,idih)
        if (NSolVariable(isol).gt.0 .and. DifferentLambdas) then
          angf = SolDihedralSize(NSolForwards(isol),idih)
          angb = SolDihedralSize(NSolBackwards(isol),idih)
        else
          angf = ang
          angb = ang
        endif
      else 
        write(printstring,*) "Unknown dihedralEnergy system type requested in calcDihedral()",ityp
        iprot = 0
        isol = 0
        ang = ZERO
        angf = ang
        angb = ang
        iparam = 0
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is probably a programming bug!")
        return
      endif

      if (iparam.eq.0) then
c       this is a null dihedral - return 0.0
        dihedralEnergy = ZERO
        dihedralEnergyF = ZERO
        dihedralEnergyB = ZERO
        return
      else if (iparam.gt.0) then
        
        dihedralEnergy = dihEnergy(iparam,ang)
        
c       if all values of lambda are the same then forwards and backwards energies
c       are same as calculated energy
        if (DifferentLambdas) then
c         if this is a variable geometry dihedral then recalculate energy for different
c         geometries
          if (angf.ne.ang) then
            dihedralEnergyF = dihEnergy(iparam,angf)
          else
            dihedralEnergyF = dihedralEnergy
          endif
              
          if (angb.ne.ang) then
            dihedralEnergyB = dihEnergy(iparam,angb)
          else
            dihedralEnergyB = dihedralEnergy
          endif
        else
          dihedralEnergyF = dihedralEnergy
          dihedralEnergyB = dihedralEnergy
        endif
        
        return
      else        
c       this is a variable dihedral - we calculate the dihedral terms
c       here as storing them all would be complicated... (and there
c       are very few variable dihedrals in a simulation)      
        iparb = VariableDihParam(-iparam,1)
        iparf = VariableDihParam(-iparam,2)
        
c       the functional forms of the dihedral energy could be very different
c       (6 terms for one, then 1 term for another!), so it is best to just
c       calculate the energies of the dihedral at lambda=0 and lambda=1 and to
c       then linearly scale the results
        
        nrg0 = dihEnergy(iparb,ang)
        nrg1 = dihEnergy(iparf,ang)        
        dihedralEnergy = perturbForceField(nrg0,nrg1,lambda)
        
c       what about forwards and backwards?
        if (DifferentLambdas) then
          nrg0 = dihEnergy(iparb,angb)
          nrg1 = dihEnergy(iparf,angb)
          dihedralEnergyB = perturbForceField(nrg0,nrg1,lambdaB)
          
          nrg0 = dihEnergy(iparb,angf)
          nrg1 = dihEnergy(iparf,angf)
          dihedralEnergyF = perturbForceField(nrg0,nrg1,lambdaF)
        else
          dihedralEnergyB = dihedralEnergy
          dihedralEnergyF = dihedralEnergy
        endif
                    
      endif
      
      end
