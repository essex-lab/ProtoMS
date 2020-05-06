      subroutine resetgbenergy
      include 'dimensions.inc'
      include 'printer.inc'
      include 'flags.inc'
      include 'enums.inc'
c################################################
c
c     Simple subroutine to recalculate the gb energies
c     correctly. Note that this is done by doing 
c     the other forcefield terms as well
c     
c     (c) Julien Michel May 2005
c################################################


      write(printstring,*) 'Calculating the total energy from scratch '
      call printLine(INFO,printstring)
      GotInitialTotalEnergy = .false.
      call checkGotTotal
      
      return
      end
