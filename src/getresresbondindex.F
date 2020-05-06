      
      integer function getResResBondIndex(iprot,sres1,sres2)
      include 'dimensions.inc'
      include 'proteins.inc'
c##################################################################
c
c     This returns the index to the array of bonded residue-residue
c     interactions - this is used by the protein non-bonded code!
c
c     (C) Christopher Woods   21 March 2002
c
c##################################################################

      integer iprot,sres1,sres2,ires1,ires2,nbnd,i
      
      if (sres1.gt.sres2) then
        ires1 = sres2
        ires2 = sres1
      else
        ires1 = sres1
        ires2 = sres2
      endif
      
      getResResBondIndex = 0
      if (ires1.eq.ires2) return
      if (.not.ResResBond(iprot,ires1,ires2)) return
      
      nbnd = NResResBonds(iprot)
      if (nbnd.le.0) return
      do i=1,nbnd
        if (ResResBondIndex(iprot,i,1).eq.ires1) then
          if (ResResBondIndex(iprot,i,2).eq.ires2) then
            getResResBondIndex = i
            return
          endif
        endif
      enddo

      return
      end
