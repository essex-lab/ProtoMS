
      double precision function DihedralHarmonicRestraintEnergy(savenrg,inrgb,inrgf)
      include 'dimensions.inc'
      include 'movelist.inc'
      include 'printer.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'extraenergies.inc'
c###############################################################################
c
c     This function calculates the total harmonic restraint energy
c     that has been applied on dihedrals
c
c     Julien Michel February 2005
c
c###############################################################################

      integer savenrg
      integer i
      double precision inrgb,inrgf
      double precision dist,distb,distf
      logical hasAtomMoved
      double precision dih,dihf,dihb
      double precision diheq,dihconst
      double precision at1c(3),at1cb(3),at1cf(3)
      double precision at2c(3),at2cb(3),at2cf(3)
      double precision at3c(3),at3cb(3),at3cf(3)
      double precision at4c(3),at4cb(3),at4cf(3)
      double precision calcGeneralDihedral
      double precision getZmatDihedral
      double precision rad2deg

      DihedralharmonicRestraintEnergy = ZERO
      inrgb = ZERO
      inrgf = ZERO
      
      if (nDihedralHarmonicRestraints.le.0) return
      
      if (allmoved) then
c       everything has moved - evaluate all of the harmonic energies
        do i=1,nDihedralHarmonicRestraints
          dih = 0.0d0
          dihf = 0.0d0
          dihb = 0.0d0
          if (zmatDihedralHarm(i)) then
            dih = getZmatDihedral(DihedralHarmAtom1(i,1),
     .            DihedralHarmAtom1(i,2),
     .            DihedralHarmAtom1(i,3),
     .            DihedralHarmAtom1(i,4),
     .            dihf,
     .            dihb)
          else
c           We need to calculate ourself the dihedral
c           First we get all the coordinates we need
c           coordinates of atX in atXc,atXcb,atXcf
            call getCoordinates(DihedralHarmAtom1(i,1),DihedralHarmAtom1(i,2),DihedralHarmAtom1(i,4),at1c,at1cb,at1cf)
            call getCoordinates(DihedralHarmAtom2(i,1),DihedralHarmAtom2(i,2),
     .        DihedralHarmAtom2(i,4),at2c,at2cb,at2cf)
            call getCoordinates(DihedralHarmAtom3(i,1),DihedralHarmAtom3(i,2),
     .        DihedralHarmAtom3(i,4),at3c,at3cb,at3cf)
            call getCoordinates(DihedralHarmAtom4(i,1),DihedralHarmAtom4(i,2),
     .        DihedralHarmAtom4(i,4),at4c,at4cb,at4cf)
c           Now calculate the dihedral
            dih = calcGeneralDihedral(at1c,at2c,at3c,at4c)
            dihf = calcGeneralDihedral(at1cf,at2cf,at3cf,at4cf)
            dihb = calcGeneralDihedral(at1cb,at2cb,at3cb,at4cb)
          endif
c          write(printstring,*) 'dih is ',rad2deg(dih),' f ',rad2deg(dihf),' b ',rad2deg(dihb)
c          call printLine(DEBUG,printstring)
          dihconst = DihedralHarmConst(i)
          diheq = DihedralHarmEq(i)
          dist = (dih - diheq)**2
          distf = (dihf - diheq)**2
          distb = (dihb - diheq)**2
          
          DihedralHarmonicRestraintEnergy = 
     .      DihedralHarmonicRestraintEnergy + dihconst*dist
          inrgf = inrgf + dihconst*distf
          inrgb = inrgb + dihconst*distb
c          write(printstring,*) 'DIH NRG ',DihedralHarmonicRestraintEnergy
c          call printLine(DEBUG,printstring)
        enddo
      else
c       loop through all of the restraints and see if they need updating
        do i=1,nDihedralHarmonicRestraints
c         If any atom has moved, calculate torsion
          if (hasAtomMoved(DihedralHarmAtom1(i,1),DihedralHarmAtom1(i,2),DihedralHarmAtom1(i,4)) 
     .      .or. hasAtomMoved(DihedralHarmAtom2(i,1),DihedralHarmAtom2(i,2),DihedralHarmAtom2(i,4)) 
     .      .or. hasAtomMoved(DihedralHarmAtom3(i,1),DihedralHarmAtom3(i,2),DihedralHarmAtom3(i,4)) 
     .      .or. hasAtomMoved(DihedralHarmAtom4(i,1),DihedralHarmAtom4(i,2),DihedralHarmAtom4(i,4)) ) 
     .      then
c            write(printstring,*) 'Some atoms may have moved recalculate'
c            call printLine(DEBUG,printstring)
            if (zmatDihedralHarm(i)) then
              dih = getZmatDihedral(DihedralHarmAtom1(i,1),
     .              DihedralHarmAtom1(i,2),
     .              DihedralHarmAtom1(i,3),
     .              DihedralHarmAtom1(i,4),
     .              dihf,
     .              dihb)
            else
              call getCoordinates(DihedralHarmAtom1(i,1),DihedralHarmAtom1(i,2),DihedralHarmAtom1(i,4),at1c,at1cb,at1cf)
              call getCoordinates(DihedralHarmAtom2(i,1),DihedralHarmAtom2(i,2),
     .          DihedralHarmAtom2(i,4),at2c,at2cb,at2cf)
              call getCoordinates(DihedralHarmAtom3(i,1),DihedralHarmAtom3(i,2),
     .          DihedralHarmAtom3(i,4),at3c,at3cb,at3cf)
              call getCoordinates(DihedralHarmAtom4(i,1),DihedralHarmAtom4(i,2),
     .          DihedralHarmAtom4(i,4),at4c,at4cb,at4cf)
c             Now calculate the dihedral
              dih = calcGeneralDihedral(at1c,at2c,at3c,at4c)
              dihf = calcGeneralDihedral(at1cf,at2cf,at3cf,at4cf)
              dihb = calcGeneralDihedral(at1cb,at2cb,at3cb,at4cb)              
            endif
            dihconst = DihedralHarmConst(i)
            diheq = DihedralHarmEq(i)
            dist = (dih - diheq)**2
            distf = (dihf - diheq)**2
            distb = (dihb - diheq)**2
            write(printstring,*) 'dih is ',rad2deg(dih),' f ',rad2deg(dihf),' b ',rad2deg(dihb), ' eq ',diheq
            call printLine(DEBUG,printstring)
            
            DihedralHarmonicRestraintEnergy = 
     .        DihedralHarmonicRestraintEnergy + dihconst*dist
            inrgf = inrgf + dihconst*distf
            inrgb = inrgb + dihconst*distb 
c            write(printstring,*) 'DIH NRG ',DihedralHarmonicRestraintEnergy
c            call printLine(DEBUG,printstring)           
          endif
        enddo
      endif      

      if (savenrg.eq.OLD) then
        oldDihHarmEnergyPart = DihedralharmonicRestraintEnergy
        oldDihHarmEnergyPartB = inrgb
        oldDihHarmEnergyPartF = inrgf
      else if (savenrg.eq.NEW) then
        newDihHarmEnergyPart = DihedralharmonicRestraintEnergy
        newDihHarmEnergyPartB = inrgb
        newDihHarmEnergyPartF = inrgf
      endif      

      return
      end
