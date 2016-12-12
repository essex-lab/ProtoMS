
      double precision function harmonicRestraintEnergy(savenrg,inrgb,inrgf)
      include 'dimensions.inc'
      include 'movelist.inc'
      include 'printer.inc'
      include 'constants.inc'
      include 'enums.inc'
      include 'extraenergies.inc'
      include 'simparams.inc'
c###############################################################################
c
c     This function calculates the total harmonic restraint energy
c     for all of the atoms that have moved.
c
c     (C) Christopher Woods, November 2004
c
c###############################################################################

      integer savenrg
      integer i
      double precision inrgb,inrgf
      double precision ic(3),icb(3),icf(3),dist,distb,distf
      logical hasAtomMoved
      
      harmonicRestraintEnergy = ZERO
      inrgb = ZERO
      inrgf = ZERO
      
      if (nHarmonicRestraints.le.0) return
      
      if (allmoved) then
c       everything has moved - evaluate all of the harmonic energies
        do i=1,nHarmonicRestraints
c         get the coordinates of the restrained atom
          call getCoordinates(HarmAtom(i,1),HarmAtom(i,2),HarmAtom(i,3),ic,icb,icf)
c         calculate the distance to the restrained point
          dist = (ic(1)-HarmPoint(i,1))**2 + (ic(2)-HarmPoint(i,2))**2 + (ic(3)-HarmPoint(i,3))**2
          distb = (icb(1)-HarmPoint(i,1))**2 + (icb(2)-HarmPoint(i,2))**2 + (icb(3)-HarmPoint(i,3))**2
          distf = (icf(1)-HarmPoint(i,1))**2 + (icf(2)-HarmPoint(i,2))**2 + (icf(3)-HarmPoint(i,3))**2
          

          if ( HarmScale(i) ) then
             harmonicRestraintEnergy = harmonicRestraintEnergy + HarmConst(i)*dist*Lambda
             inrgb = inrgb + HarmConst(i)*distb*LambdaB
             inrgf = inrgf + HarmConst(i)*distf*LambdaF
          else
             harmonicRestraintEnergy = harmonicRestraintEnergy + HarmConst(i)*dist
             inrgb = inrgb + HarmConst(i)*distb
             inrgf = inrgf + HarmConst(i)*distf
          endif 
        enddo
      else
c       loop through all of the restraints and see if they need updating
        do i=1,nHarmonicRestraints
          if (hasAtomMoved(HarmAtom(i,1),HarmAtom(i,2),HarmAtom(i,3))) then
            call getCoordinates(HarmAtom(i,1),HarmAtom(i,2),HarmAtom(i,3),ic,icb,icf)

c           calculate the distance to the restrained point
            dist = (ic(1)-HarmPoint(i,1))**2 + (ic(2)-HarmPoint(i,2))**2 + (ic(3)-HarmPoint(i,3))**2
            distb = (icb(1)-HarmPoint(i,1))**2 + (icb(2)-HarmPoint(i,2))**2 + (icb(3)-HarmPoint(i,3))**2
            distf = (icf(1)-HarmPoint(i,1))**2 + (icf(2)-HarmPoint(i,2))**2 + (icf(3)-HarmPoint(i,3))**2
          

            if ( HarmScale(i) ) then
               harmonicRestraintEnergy = harmonicRestraintEnergy + HarmConst(i)*dist*Lambda
               inrgb = inrgb + HarmConst(i)*distb*LambdaB
               inrgf = inrgf + HarmConst(i)*distf*LambdaF
            else
               harmonicRestraintEnergy = harmonicRestraintEnergy + HarmConst(i)*dist
               inrgb = inrgb + HarmConst(i)*distb
               inrgf = inrgf + HarmConst(i)*distf
            endif 
          endif
        enddo
      endif      
        
      if (savenrg.eq.OLD) then
        oldHarmEnergyPart = harmonicRestraintEnergy
        oldHarmEnergyPartB = inrgb
        oldHarmEnergyPartF = inrgf
        write(printstring,*) 'old harmNrg ?',oldHarmEnergyPart,' F ',oldHarmEnergyPartF,
     .    ' B ',oldHarmEnergyPartB
        call printLine(DEBUG,printstring)      
      else if (savenrg.eq.NEW) then
        newHarmEnergyPart = harmonicRestraintEnergy
        newHarmEnergyPartB = inrgb
        newHarmEnergyPartF = inrgf
        write(printstring,*) 'new harmNrg ?',newHarmEnergyPart,' F ',newHarmEnergyPartF,
     .    ' B ',newHarmEnergyPartB
        call printLine(DEBUG,printstring)      
      endif      

      return
      end
