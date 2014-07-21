      subroutine switchBornRadii(confA,confB)
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
c     This subroutine put the Born radii data of confA 
c     into '...' and confB into 'Old...'
c     
c     (c) Julien Michel May 2005
c###########################################################
      integer i,j,confA,confB

c      write(printstring,*) 'STORING COORDINATES of conf ',conf
c      call printLine(DEBUG,printstring)

      if (NProteins .gt. 0) then
        do i=1,NProteins
          do j=1,NZAtoms(i)
            ZProtContBorn(i,j) = StoreZProtContBorn(i,j,confA) 
            ZProtContBornF(i,j) = StoreZProtContBornF(i,j,confA) 
            ZProtContBornB(i,j) = StoreZProtContBornB(i,j,confA) 
            OldZProtContBorn(i,j) = StoreOldZProtContBorn(i,j,confB) 
            OldZProtContBornF(i,j) = StoreOldZProtContBornF(i,j,confB) 
            OldZProtContBornB(i,j) = StoreOldZProtContBornB(i,j,confB) 
            OldZProtEffBorn(i,j) = StoreOldZProtEffBorn(i,j,confB) 
            OldZProtEffBornF(i,j) = StoreOldZProtEffBornF(i,j,confB) 
            OldZProtEffBornB(i,j) = StoreOldZProtEffBornB(i,j,confB) 
          enddo
        enddo
      endif

      if (NSolutes.gt.0) then
        do i=1,NSolutes
          do j=1,NSolAtoms(i)
            ZSolContBorn(i,j) = StoreZSolContBorn(i,j,confA) 
            ZSolContBornF(i,j) = StoreZSolContBornF(i,j,confA) 
            ZSolContBornB(i,j) = StoreZSolContBornB(i,j,confA) 
            OldZSolContBorn(i,j) = StoreOldZSolContBorn(i,j,confB) 
            OldZSolContBornF(i,j) = StoreOldZSolContBornF(i,j,confB) 
            OldZSolContBornB(i,j) = StoreOldZSolContBornB(i,j,confB) 
            OldZSolEffBorn(i,j) = StoreOldZSolEffBorn(i,j,confB)
            OldZSolEffBornF(i,j) = StoreOldZSolEffBornF(i,j,confB)
            OldZSolEffBornB(i,j) = StoreOldZSolEffBornB(i,j,confB)
          enddo
        enddo
      endif

      if (NSolvents .gt. 0) then
        do i=1,NSolvents
          do j=1,NSvnAtoms(i)
            ZSvnContBorn(i,j) = StoreZSvnContBorn(i,j,confA) 
            ZSvnContBornF(i,j) = StoreZSvnContBornF(i,j,confA) 
            ZSvnContBornB(i,j) = StoreZSvnContBornB(i,j,confA)
            OldZSvnContBorn(i,j) = StoreOldZSvnContBorn(i,j,confB) 
            OldZSvnContBornF(i,j) = StoreOldZSvnContBornF(i,j,confB) 
            OldZSvnContBornB(i,j) = StoreOldZSvnContBornB(i,j,confB)
            OldZSvnEffBorn(i,j) = StoreOldZSvnEffBorn(i,j,confB) 
            OldZSvnEffBornF(i,j) = StoreOldZSvnEffBornF(i,j,confB) 
            OldZSvnEffBornB(i,j) = StoreOldZSvnEffBornB(i,j,confB)
          enddo
        enddo
      endif

      return
      end
