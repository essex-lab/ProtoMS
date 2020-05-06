      subroutine StoreBornData(pot)
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
c     This subroutine save the Born radii of the entire 
c     system
c
c     (c) Julien Michel May 2005
c###########################################################
      integer i,j,pot

c      write(printstring,*) 'STORING COORDINATES of conf ',conf
c      call printLine(DEBUG,printstring)

      if (NProteins .gt. 0) then
        do i=1,NProteins
          do j=1,NZAtoms(i)
            StoreZProtEffBorn(i,j,pot) = ZProtEffBorn(i,j)
            StoreZProtEffBornF(i,j,pot) = ZProtEffBornF(i,j)
            StoreZProtEffBornB(i,j,pot) = ZProtEffBornB(i,j)
            StoreZProtContBorn(i,j,pot) = ZProtContBorn(i,j)
            StoreZProtContBornF(i,j,pot) = ZProtContBornF(i,j)
            StoreZProtContBornB(i,j,pot) = ZProtContBornB(i,j)
            StoreOldZProtContBorn(i,j,pot) = OldZProtContBorn(i,j)
            StoreOldZProtContBornF(i,j,pot) = OldZProtContBornF(i,j)
            StoreOldZProtContBornB(i,j,pot) = OldZProtContBornB(i,j)
            StoreOldZProtEffBorn(i,j,pot) = OldZProtEffBorn(i,j)
            StoreOldZProtEffBornF(i,j,pot) = OldZProtEffBornF(i,j)
            StoreOldZProtEffBornB(i,j,pot) = OldZProtEffBornB(i,j)
          enddo
        enddo
      endif

      if (NSolutes.gt.0) then
        do i=1,NSolutes
          do j=1,NSolAtoms(i)
            StoreZSolEffBorn(i,j,pot) = ZSolEffBorn(i,j)
            StoreZSolEffBornF(i,j,pot) = ZSolEffBornF(i,j)
            StoreZSolEffBornB(i,j,pot) = ZSolEffBornB(i,j)
            StoreZSolContBorn(i,j,pot) = ZSolContBorn(i,j)
            StoreZSolContBornF(i,j,pot) = ZSolContBornF(i,j)
            StoreZSolContBornB(i,j,pot) = ZSolContBornB(i,j)
            StoreOldZSolContBorn(i,j,pot) = OldZSolContBorn(i,j)
            StoreOldZSolContBornF(i,j,pot) = OldZSolContBornF(i,j)
            StoreOldZSolContBornB(i,j,pot) = OldZSolContBornB(i,j)
            StoreOldZSolEffBorn(i,j,pot) = OldZSolEffBorn(i,j)
            StoreOldZSolEffBornF(i,j,pot) = OldZSolEffBornF(i,j)
            StoreOldZSolEffBornB(i,j,pot) = OldZSolEffBornB(i,j)
          enddo
        enddo
      endif

      if (NSolvents .gt. 0) then
        do i=1,NSolvents
          do j=1,NSvnAtoms(i)
            StoreZSvnEffBorn(i,j,pot) = ZSvnEffBorn(i,j)
            StoreZSvnEffBornF(i,j,pot) = ZSvnEffBornF(i,j)
            StoreZSvnEffBornB(i,j,pot) = ZSvnEffBornB(i,j)
            StoreZSvnContBorn(i,j,pot) = ZSvnContBorn(i,j)
            StoreZSvnContBornF(i,j,pot) = ZSvnContBornF(i,j)
            StoreZSvnContBornB(i,j,pot) = ZSvnContBornB(i,j)
            StoreOldZSvnContBorn(i,j,pot) = OldZSvnContBorn(i,j)
            StoreOldZSvnContBornF(i,j,pot) = OldZSvnContBornF(i,j)
            StoreOldZSvnContBornB(i,j,pot) = OldZSvnContBornB(i,j)
            StoreOldZSvnEffBorn(i,j,pot) = OldZSvnEffBorn(i,j)
            StoreOldZSvnEffBornF(i,j,pot) = OldZSvnEffBornF(i,j)
            StoreOldZSvnEffBornB(i,j,pot) = OldZSvnEffBornB(i,j)
          enddo
        enddo
      endif

      return
      end
