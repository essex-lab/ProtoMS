
      double precision function hardwallEnergy(savenrg,inrgb,inrgf)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'constants.inc'
      include 'enums.inc'
      include 'extraenergies.inc'
      include 'solvents.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      
c###############################################################################
c
c     This calculates the hardwall potential applied to each atom
c     If the atoms are within a certain radius defined by HardwallRadius (previously
c        defined by irad) then they will experience a large potential
c     This applies to all but the molecule of interest, which will experience a large
c        potential if it goes outside of the defined radius
c
c     J. Taylor July 2005
c
c     Based on the ProtoMS version created by C. Woods and editted by C. Barillari on 05-05-2004  
c
c###############################################################################

      integer savenrg
      integer i,j,k,jj
      double precision inrgb,inrgf
      double precision RadiusSq,idist,idist2
     
      
      hardwallEnergy = ZERO
      inrgb = ZERO
      inrgf = ZERO
      
      if (nHardwallPotential.le.0) return

      do i=1,nHardwallPotential
c     The hard-wall potential is centered on a point (HardwallCenter) which is given as input. It corresponds to the
c     original cartesian coordinates of the oxygen of the water. The radius on which the HW potential acts around the
c     defined point is user-defined (HardwallRadius). 
c     
c     The distance between the solvent atoms and the defined point is calculated and if the distance is
c     smaller than HardwallRadius the energy is very high (HardwallEnergy)
        
      RadiusSq = HardwallRadius(i)**2
          
      if (NSolvents.gt.0) then
        do 200 j=1,NSolvents
      do 201 k=1,NSvnAtoms(j) 
        idist  = ((HardwallCenter(i,1)-ZSvnCoords(j,k,1))**2 + (HardwallCenter(i,2)-ZSvnCoords(j,k,2))**2 
     .                + (HardwallCenter(i,3)-ZSvnCoords(j,k,3))**2)
            if (idist .le. RadiusSq) then
                hardwallEnergy = hardwallEnergy + HardwallextraEnergy(i)
                inrgb = hardwallEnergy
                inrgf = hardwallEnergy
c            call printLine(DEBUG,"Hardwall applied to some solvent atom, MC fails") 
            goto 300
            endif 
201      continue
200     continue             
      endif


c     Then the distance between each atom of the protein and the reference point is calculated and compared to the
c     user defined radius.  
          
        if (NProteins.gt.0) then
        do 210 j=1,NProteins
        do 220 k=1,NZAtoms(j)
        idist = ((HardwallCenter(i,1)-ZCoords(j,k,1))**2 + (HardwallCenter(i,2)-ZCoords(j,k,2))**2 
     .            + (HardwallCenter(i,3)-ZCoords(j,k,3))**2)
            if (idist .le. RadiusSq) then
                hardwallEnergy = hardwallEnergy + HardwallextraEnergy(i)
                inrgb = hardwallEnergy
                inrgf = hardwallEnergy
c            call printLine(RESULTS,"Hardwall applied to some protein atom, MC fails")
              goto 300
            endif
220     continue
210     continue  
        endif

c    This should work for any kind of solute
c    Also, only the given atom of the solute should feel the hardwall within the residue

c    The solute is more complicated. The HW is built around a water molecule, and has the double function of
c    keeping the water inside and the rest outside.
c    Differentiate between the the solutes using the id assigned by processhardwallline.F and addid.F 
c    When the solute is anything but the molecules of interest the energy must be high if the
c    distance between the its atoms and the reference point is less than the user-defined distance.
c    When the solute is the water of interest, the opposite must be verified.     
               
c        write(6,*) "HardwallMolecule",Hardwallmolecule(i)
c   	write(6,*) "HardwallCenter",HardwallCenter(i,1),HardwallCenter(i,2),HardwallCenter(i,3)

        if (NSolutes.gt.0) then
          do 235 j=1,NSolutes
c           for the molecule of interest, we check whether its geometry center is in the region of the hardwall
            if (j.eq.HardwallMolecule(i)) then
              idist = ((HardwallCenter(i,1)-SolGeomCentre(j,1))**2 + (HardwallCenter(i,2)-SolGeomCentre(j,2))**2 
     .              + (HardwallCenter(i,3)-SolGeomCentre(j,3))**2)
c              write(6,*) "SolGeomCentre",SolGeomCentre(j,1),SolGeomCentre(j,2),SolGeomCentre(j,3)
c              write(6,*) "Distance",(idist)**0.5
              if (idist.ge.RadiusSq) then
                hardwallEnergy  = hardwallEnergy + HardwallextraEnergy(i)
                inrgb = hardwallEnergy
                inrgf = hardwallEnergy
                goto 300
              endif
            else
              do 245 k=1,NSolAtoms(j)
                if (.not.IgnoreSolute(j)) then
                  idist = ((HardwallCenter(i,1)-ZSolCoords(j,k,1))**2 + (HardwallCenter(i,2)-ZSolCoords(j,k,2))**2 
     .              + (HardwallCenter(i,3)-ZSolCoords(j,k,3))**2)
                else
                  if (idist.le.RadiusSq) then
                    hardwallEnergy  = hardwallEnergy + HardwallextraEnergy(i)
                    inrgb = hardwallEnergy
                    inrgf = hardwallEnergy
                    goto 300
                  endif
                endif  
245           continue          
            endif  
235       continue
c245       continue
c235       continue
        endif 

      enddo
   
       
300   if (savenrg.eq.OLD) then
        oldHardwallEnergyPart = hardwallEnergy
        oldHardwallEnergyPartB = inrgb
        oldHardwallEnergyPartF = inrgf
c       write(printstring,*) 'old hardwallNrg ?',oldHardwallEnergyPart,' F ',oldHardwallEnergyPartF,' B ',oldHardwallEnergyPartB
c       call printLine(DEBUG,printstring)      
      else if (savenrg.eq.NEW) then
        newHardwallEnergyPart = hardwallEnergy
        newHardwallEnergyPartB = inrgb
        newHardwallEnergyPartF = inrgf
c       write(printstring,*) 'new hardwallNrg ?',newHardwallEnergyPart,' F ',newHardwallEnergyPartF,' B ',newHardwallEnergyPartB
c       call printLine(DEBUG,printstring)      
      endif      


      return
      end
 
 
