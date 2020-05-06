      subroutine FlagMovedBits
      include 'dimensions.inc'
      include 'movelist.inc'
      include 'solutes.inc'
      include 'quickpot.inc'
c##############################################
c     Add every elements of splitMovedList
c     onto MovedList
c
c
c     (c) Julien Michel April 2005
c
c##############################################
      integer i,isol,iprot,ires,isvn

c     For the solutes
      do i=1,nSolSplitMoved
        isol = SplitSolMovedList(i)
        call setSoluteMoved(isol)
      enddo
c     For the residues
      do i=1,nResSplitMoved
        iprot = SplitResMovedList(i,1)
        ires = SplitResMovedList(i,2)
        call setResidueMoved(iprot,ires)
      enddo
c     For the solvents
      do i=1,nSvnSplitMoved
        isvn = SplitSvnMovedList(i)
        call setSolventMoved(isvn)
      enddo
      return
      end
