      subroutine fakeSimulate()
      include 'enums.inc'
c##############################################################
c
c     This performs one step of simulation, without doing
c     anything other than adding the energies to the 
c     averages. This is useful when for a developer who
c     has to debug energy components (like me!)
c
c     (C) Christopher Woods, October 2004
c
c##############################################################

      double precision inrg,inrgb,inrgf
      double precision totalEnergy,extraEnergy
      
c     need to set all moved so that the total extra energy is calculated correctly
      call setAllMoved     

      inrg = totalEnergy(OLD,inrgb,inrgf)+extraEnergy(OLD,inrgb,inrgf)
c     can now clear the move list
      call clearMoveList  

      call calcAverageEnergy
      return
      end
