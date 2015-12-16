      
      subroutine singlePointEnergy
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
c###########################################################
c
c     This routine calculates and prints the energy of a single
c     conformation of the system
c
c     (C) Christopher Woods, November 2004
c
c###########################################################

      double precision tempval,nrgb,nrgf,inrgb,inrgf
      double precision totalEnergy,extraEnergy
      double precision tstart,tend,telapsed

c     need to set all moved so that the total extra energy is calculated correctly
      call setAllMoved      
      
      call cpu_time(tstart)

      tempval = totalEnergy(DONTCARE,nrgb,nrgf) + extraEnergy(DONTCARE,inrgb,inrgf)
 
      call cpu_time(tend)
      telapsed = tend - tstart
 
      nrgb = nrgb + inrgb
      nrgf = nrgf + inrgf

c     can now clear the move list
      call clearMoveList      
            
      write(printstring,*) "Singlepoint energy = ",tempval,' fwd ',nrgf,' bwd ',nrgb
      call printLine(INFO,printstring)
      call printLine(SPENERGY,printstring)

      write(printstring,"('Energy calculation took ',EN8.1E1,' s')") telapsed
      call printLine(INFO,printstring)

      return
      end
