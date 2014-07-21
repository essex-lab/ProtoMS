      subroutine moveSideChain(iprot,ires)
      include 'dimensions.inc'
      include 'proteins.inc'
c################################################################
c
c     This routine moves the sidechain of residue
c     ires in protein iprot
c
c     (C) Christopher Woods    6 February 2002
c
c################################################################

      integer iprot,ires,i,j,nsc
      double precision ranval,ranRange

c     move the side chain atoms...
      nsc = NSCAtoms(iprot,ires)
      if (nsc.gt.0) then
        do j=1,nsc
          do i=1,3
            if (ZMovSC(iprot,ires,j,i)) then
              ranval = ranRange(-ZDeltaSC(iprot,ires,j,i),ZDeltaSC(iprot,ires,j,i))
              ZValSC(iprot,ires,j,i) = ZValSC(iprot,ires,j,i) + ranval
            endif
          enddo
        enddo
      endif

      return
      end
