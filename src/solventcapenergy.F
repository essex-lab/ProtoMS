      
      double precision function solventCapEnergy(savenrg,inrgb,inrgf)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'enums.inc'
      include 'movelist.inc'
      include 'flags.inc'
      include 'constants.inc'
      include 'extraenergies.inc'
      include 'simparams.inc'
      include 'printer.inc'
c################################################################
c
c     This calculates the energy of the solvent cap for
c     the parts of the system that have moved.
c
c     (C) Christopher Woods, November 2004
c
c#################################################################

      integer savenrg,i,isvn
      double precision inrgb,inrgf
      double precision idist,caprad2
            
      solventCapEnergy = ZERO
      inrgb = ZERO
      inrgf = ZERO
      
      if (NSolvents.le.0 .or. .not.SolventCap) return
      
      caprad2 = capradius**2
      
      if (allmoved) then
c       calculate the solvent cap energy of all solvent molecules
        do isvn=1,NSolvents
          idist = (SvnCentre(isvn,1)-cap(1))**2 + (SvnCentre(isvn,2)-cap(2))**2
     .            + (SvnCentre(isvn,3)-cap(3))**2
     
          if (idist.gt.caprad2) then
            solventCapEnergy = solventCapEnergy + 
     .                  capconstant*( dsqrt(idist) - capradius)**2
          endif
        enddo
        
c       cap energy is constant with respect to lambda
        inrgf = solventCapEnergy
        inrgb = solventCapEnergy        

      else
c       loop over all the moved solvent molecules and calculate their
c       cap energy
        if (nSvnMoved.gt.0) then
          do i=1,nSvnMoved
            isvn = svnMovedList(i)

            idist = (SvnCentre(isvn,1)-cap(1))**2 + (SvnCentre(isvn,2)-cap(2))**2
     .              + (SvnCentre(isvn,3)-cap(3))**2
     
            if (idist.gt.caprad2) then
              solventCapEnergy = solventCapEnergy + 
     .                    capconstant*( dsqrt(idist) -capradius)**2
c            write(printstring,*) 'isvn ',isvn,' idist ',idist,' caprad2 ',caprad2
c            call printLine(DEBUG,printstring)
c            write(printstring,*) 'capconst ',capconstant,' caprad ',capradius
c            call printLinE(DEBUG,printstring)
c            write(printstring,*) 'solventCapEnergy',solventCapEnergy
c            call printLine(DEBUG,printstring)
c            stop
            endif

          enddo
c         cap energy is constant with respect to lambda
          inrgf = solventCapEnergy
          inrgb = solventCapEnergy          
        
        endif  
      endif
c     see if we need to save the energies
      if (savenrg.eq.OLD) then
        oldSvnCapEnergyPart = solventCapEnergy
c        write(printstring,*) 'oldsolventcapNRG ',solventcapEnergy
c        call printLine(DEBUG,printstring)
      else if (savenrg.eq.NEW) then
        newSvnCapEnergyPart = solventCapEnergy
c        write(printstring,*) 'newsolventcapNRG ',solventcapEnergy
c        call printLine(DEBUG,printstring)
      endif                

      return
      end
