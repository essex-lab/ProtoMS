      subroutine RevertBornData(pot)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'solvents.inc'
      include 'flags.inc'
      include 'energies.inc'
      include 'printer.inc'
      include 'gb.inc'
      include 'quickpot.inc'
c############################################################
c     This subroutine reverts the Born radii of the entire 
c     system for a given conf
c
c     (c) Julien Michel May 2005
c###########################################################
      integer i,j,pot

c      write(printstring,*) 'STORING COORDINATES of pot ',pot
c      call printLine(DEBUG,printstring)

      if (NProteins .gt. 0) then
        do i=1,NProteins
          do j=1,NZAtoms(i)
            ZProtContBorn(i,j) = StoreZProtContBorn(i,j,pot) 
            ZProtContBornF(i,j) = StoreZProtContBornF(i,j,pot) 
            ZProtContBornB(i,j) = StoreZProtContBornB(i,j,pot) 
            ZProtEffBorn(i,j) = StoreZProtEffBorn(i,j,pot) 
            ZProtEffBornF(i,j) = StoreZProtEffBornF(i,j,pot) 
            ZProtEffBornB(i,j) = StoreZProtEffBornB(i,j,pot) 
            OldZProtContBorn(i,j) = StoreOldZProtContBorn(i,j,pot) 
            OldZProtContBornF(i,j) = StoreOldZProtContBornF(i,j,pot) 
            OldZProtContBornB(i,j) = StoreOldZProtContBornB(i,j,pot) 
            OldZProtEffBorn(i,j) = StoreOldZProtEffBorn(i,j,pot) 
            OldZProtEffBornF(i,j) = StoreOldZProtEffBornF(i,j,pot) 
            OldZProtEffBornB(i,j) = StoreOldZProtEffBornB(i,j,pot) 
          enddo
        enddo
      endif

      if (NSolutes.gt.0) then
        do i=1,NSolutes
          do j=1,NSolAtoms(i)
            ZSolContBorn(i,j) = StoreZSolContBorn(i,j,pot) 
            ZSolContBornF(i,j) = StoreZSolContBornF(i,j,pot) 
            ZSolContBornB(i,j) = StoreZSolContBornB(i,j,pot) 
            ZSolEffBorn(i,j) = StoreZSolEffBorn(i,j,pot) 
            ZSolEffBornF(i,j) = StoreZSolEffBornF(i,j,pot) 
            ZSolEffBornB(i,j) = StoreZSolEffBornB(i,j,pot) 
            OldZSolContBorn(i,j) = StoreOldZSolContBorn(i,j,pot) 
            OldZSolContBornF(i,j) = StoreOldZSolContBornF(i,j,pot) 
            OldZSolContBornB(i,j) = StoreOldZSolContBornB(i,j,pot) 
            OldZSolEffBorn(i,j) = StoreOldZSolEffBorn(i,j,pot)
            OldZSolEffBornF(i,j) = StoreOldZSolEffBornF(i,j,pot)
            OldZSolEffBornB(i,j) = StoreOldZSolEffBornB(i,j,pot)
          enddo
        enddo
      endif

      if (NSolvents .gt. 0) then
        do i=1,NSolvents
          do j=1,NSvnAtoms(i)
            ZSvnContBorn(i,j) = StoreZSvnContBorn(i,j,pot) 
            ZSvnContBornF(i,j) = StoreZSvnContBornF(i,j,pot) 
            ZSvnContBornB(i,j) = StoreZSvnContBornB(i,j,pot)
            ZSvnEffBorn(i,j) = StoreZSvnEffBorn(i,j,pot) 
            ZSvnEffBornF(i,j) = StoreZSvnEffBornF(i,j,pot) 
            ZSvnEffBornB(i,j) = StoreZSvnEffBornB(i,j,pot)
            OldZSvnContBorn(i,j) = StoreOldZSvnContBorn(i,j,pot) 
            OldZSvnContBornF(i,j) = StoreOldZSvnContBornF(i,j,pot) 
            OldZSvnContBornB(i,j) = StoreOldZSvnContBornB(i,j,pot)
            OldZSvnEffBorn(i,j) = StoreOldZSvnEffBorn(i,j,pot) 
            OldZSvnEffBornF(i,j) = StoreOldZSvnEffBornF(i,j,pot) 
            OldZSvnEffBornB(i,j) = StoreOldZSvnEffBornB(i,j,pot)
          enddo
        enddo
      endif

      return
      end
