
      double precision function flatbottomRestraintEnergy(savenrg,inrgb,inrgf)
      include 'dimensions.inc'
      include 'movelist.inc'
      include 'printer.inc'
      include 'constants.inc'
      include 'enums.inc'
      include 'extraenergies.inc'
c###############################################################################
c
c     This function calculates the total flat bottom harmonic restraint energy
c     for all of the atoms that have moved.
c
c     Julien Michel February 2005 based on harmonic restraint energy
c
c###############################################################################

      integer savenrg
      integer i
      double precision inrgb,inrgf
      double precision width,const
      double precision ic(3),icb(3),icf(3),dist,distb,distf
      logical hasAtomMoved
      
      flatbottomRestraintEnergy = ZERO
      inrgb = ZERO
      inrgf = ZERO
      
      if (nFlatBottomRestraints.le.0) return
      
      if (allmoved) then
c       everything has moved - evaluate all of the flatbottom energies
        do i=1,nFlatBottomRestraints
c         get the coordinates of the restrained atom
          call getCoordinates(FBotAtom(i,1),FBotAtom(i,2),FBotAtom(i,3),ic,icb,icf)
c         calculate the distance to the restrained point
          dist = sqrt((ic(1)-FBotPoint(i,1))**2 + (ic(2)-FBotPoint(i,2))**2 + (ic(3)-FBotPoint(i,3))**2)
          distb = sqrt((icb(1)-FBotPoint(i,1))**2 + (icb(2)-FBotPoint(i,2))**2 + (icb(3)-FBotPoint(i,3))**2)
          distf = sqrt((icf(1)-FBotPoint(i,1))**2 + (icf(2)-FBotPoint(i,2))**2 + (icf(3)-FBotPoint(i,3))**2)
c         substract width of flat restraint
          width = FBotWidth(i)
          const = FBotConst(i)
c          write(printstring,*) ' dist is ',dist
c          call printLine(DEBUG,printstring)
          dist = dist - width
          distb = distb - width
          distf = distf - width
c          write(printstring,*) 'Point ',FBotPoint(i,1),FBotPoint(i,2),FBotPoint(i,3)
c          call printLine(DEBUG,printstring)
c          write(printstring,*) 'coords ',ic(1),ic(2),ic(3)
c          call printLine(DEBUG,printstring)
c          write(printstring,*) 'width ',width,' dist -witdh ',dist,' const ',FBotConst(i)
c          call printLine(DEBUG,printstring)
          if (dist .gt. ZERO) flatbottomRestraintEnergy = flatbottomRestraintEnergy + const*dist**2
          if (distf .gt. ZERO) inrgb = inrgb + const*distb**2
          if (distb .gt. ZERO) inrgf = inrgf + const*distf**2
        enddo
      else
c       loop through all of the restraints and see if they need updating
        do i=1,nFlatBottomRestraints
          if (hasAtomMoved(FBotAtom(i,1),FBotAtom(i,2),FBotAtom(i,3))) then
            call getCoordinates(FBotAtom(i,1),FBotAtom(i,2),FBotAtom(i,3),ic,icb,icf)

c           calculate the distance to the restrained point
            dist = sqrt((ic(1)-FBotPoint(i,1))**2 + (ic(2)-FBotPoint(i,2))**2 + (ic(3)-FBotPoint(i,3))**2)
            distb = sqrt((icb(1)-FBotPoint(i,1))**2 + (icb(2)-FBotPoint(i,2))**2 + (icb(3)-FBotPoint(i,3))**2)
            distf = sqrt((icf(1)-FBotPoint(i,1))**2 + (icf(2)-FBotPoint(i,2))**2 + (icf(3)-FBotPoint(i,3))**2)
c           substract width of flat restraint
            width = FBotWidth(i)
            const = FBotConst(i)
            dist = dist - width
            distb = distb - width
            distf = distf - width
            if (dist .gt. ZERO) flatbottomRestraintEnergy = flatbottomRestraintEnergy + const*dist**2
            if (distf .gt. ZERO) inrgb = inrgb + const*distb**2
            if (distb .gt. ZERO) inrgf = inrgf + const*distf**2
          endif
        enddo
      endif      
        
      if (savenrg.eq.OLD) then
        oldFBotEnergyPart = flatbottomRestraintEnergy
        oldFBotEnergyPartB = inrgb
        oldFBotEnergyPartF = inrgf
c        write(printstring,*) 'old flatbotNrg ?',oldFBotEnergyPart,' F ',oldFBotEnergyPartF,
c     .    ' B ',oldFBotEnergyPartB
c        call printLine(DEBUG,printstring)      
      else if (savenrg.eq.NEW) then
        newFBotEnergyPart = flatbottomRestraintEnergy
        newFBotEnergyPartB = inrgb
        newFBotEnergyPartF = inrgf
c        write(printstring,*) 'new fbot nrg  ?',newFBotEnergyPart,' F ',newFBotEnergyPartF,
c     .    ' B ',newFBotEnergyPartB
c        call printLine(DEBUG,printstring)      
      endif      

      return
      end
