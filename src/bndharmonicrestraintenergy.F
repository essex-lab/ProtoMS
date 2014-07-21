
      double precision function bndharmonicRestraintEnergy(savenrg,inrgb,inrgf)
      include 'dimensions.inc'
      include 'movelist.inc'
      include 'printer.inc'
      include 'constants.inc'
      include 'enums.inc'
      include 'flags.inc'
      include 'simparams.inc'
      include 'extraenergies.inc'
c###############################################################################
c
c     This function calculates the total bond harmonic restraint energy
c     for all of the atoms that have moved.
c
c     Julien Michel Jan 2007
c
c###############################################################################

      integer savenrg
      integer i
      double precision inrgb,inrgf
      double precision ic1(3),icb1(3),icf1(3),ic2(3),icb2(3),icf2(3),dist,distb,distf
      logical hasAtomMoved
      
      bndharmonicRestraintEnergy = ZERO
      inrgb = ZERO
      inrgf = ZERO
      
      if (nbndHarmonicRestraints.le.0) return
      
      if (allmoved) then
c       everything has moved - evaluate all of the harmonic energies
        do i=1,nbndHarmonicRestraints
c         get the coordinates of the first solvent 
          call getCoordinates(BndHarmAtom1(i,1),BndHarmAtom1(i,2),BndHarmAtom1(i,3),ic1,icb1,icf1)
          call getCoordinates(BndHarmAtom2(i,1),BndHarmAtom2(i,2),BndHarmAtom2(i,3),ic2,icb2,icf2)

c          calculate the distance between these two atoms
          dist = (ic1(1)-ic2(1))**2 + (ic1(2)-ic2(2))**2 + (ic1(3)-ic2(3))**2
          distb = (icb1(1)-icb2(1))**2 + (icb1(2)-icb2(2))**2 + (icb1(3)-icb2(3))**2
          distf =(icf1(1)-icf2(1))**2 + (icf1(2)-icf2(2))**2 + (icf1(3)-icf2(3))**2

          dist = dsqrt(dist)
          distb = dsqrt(distb)
          distf = dsqrt(distf)
c         Now substract the equilibrium distance

c          write(printstring,*) 'DIST spring',dist,' EQ ',BndHarmEq(i)
c          call printLine(DEBUG,printstring)

          dist =  ( dist  - BndHarmEq(i) )**2
          distb = ( distb - BndHarmEq(i) )**2
          distf = ( distf - BndHarmEq(i) )**2

        
          bndharmonicRestraintEnergy = bndharmonicRestraintEnergy + BndHarmConst(i)*dist
          inrgb = inrgb + BndHarmConst(i)*distb
          inrgf = inrgf + BndHarmConst(i)*distf
          
c          write(printstring,*) 'DIST ',dist,' nrg ',BndHarmConst(i)*dist
c          call printLine(DEBUG,printstring)

c          inrgb = inrgb + BndHarmConst(i)*distb
c          inrgf = inrgf + BndHarmConst(i)*distf
        enddo
      else
c       loop through all of the restraints and see if they need updating
        do i=1,nbndHarmonicRestraints
          if (hasAtomMoved(BndHarmAtom1(i,1),BndHarmAtom1(i,2),BndHarmAtom1(i,3)) .or. 
     .      hasAtomMoved(BndHarmAtom2(i,1),BndHarmAtom2(i,2),BndHarmAtom2(i,3))  )  then
c         get the coordinates of the first solvent 
          call getCoordinates(BndHarmAtom1(i,1),BndHarmAtom1(i,2),BndHarmAtom1(i,3),ic1,icb1,icf1)
          call getCoordinates(BndHarmAtom2(i,1),BndHarmAtom2(i,2),BndHarmAtom2(i,3),ic2,icb2,icf2)

c          calculate the distance between these two atoms
          dist = (ic1(1)-ic2(1))**2 + (ic1(2)-ic2(2))**2 + (ic1(3)-ic2(3))**2
          distb = (icb1(1)-icb2(1))**2 + (icb1(2)-icb2(2))**2 + (icb1(3)-icb2(3))**2
          distf =(icf1(1)-icf2(1))**2 + (icf1(2)-icf2(2))**2 + (icf1(3)-icf2(3))**2

          dist = dsqrt(dist)
          distb = dsqrt(distb)
          distf = dsqrt(distf)
c         Now substract the equilibrium distance
c          write(printstring,*) 'DIST spring',dist,' EQ ',BndHarmEq(i)
c          call printLine(DEBUG,printstring)

          dist =  ( dist  - BndHarmEq(i) )**2
          distb = ( distb - BndHarmEq(i) )**2
          distf = ( distf - BndHarmEq(i) )**2

c          write(printstring,*) 'DIST ',dist,' nrg ',BndHarmConst(i)*dist
c          call printLine(DEBUG,printstring)
        
          bndharmonicRestraintEnergy = bndharmonicRestraintEnergy + BndHarmConst(i)*dist
          inrgb = inrgb + BndHarmConst(i)*distb
          inrgf = inrgf + BndHarmConst(i)*distf
          endif
        enddo
      endif      
        

      if (savenrg.eq.OLD) then
        oldBndHarmEnergyPart = BndharmonicRestraintEnergy
        oldBndHarmEnergyPartB = inrgb
        oldBndHarmEnergyPartF = inrgf
c       write(printstring,*) 'old BndharmNrg ?',oldBndHarmEnergyPart,
c    .    ' F ',oldBndHarmEnergyPartF,
c    .    ' B ',oldBndHarmEnergyPartB
c       call printLine(DEBUG,printstring)      
      else if (savenrg.eq.NEW) then
        newBndHarmEnergyPart = BndharmonicRestraintEnergy
        newBndHarmEnergyPartB = inrgb
        newBndHarmEnergyPartF = inrgf
c       write(printstring,*) 'new BndharmNrg ?',newBndHarmEnergyPart,
c    .    ' F ',newBndHarmEnergyPartF,
c    .    ' B ',newBndHarmEnergyPartB
c       call printLine(DEBUG,printstring)      
      endif      

      return
      end
