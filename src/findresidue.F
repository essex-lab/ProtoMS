      
      integer function findResidue(iprot,ires)
      include 'dimensions.inc'
      include 'proteins.inc'
c############################################################
c
c     This finds and returns the index of residue 'ires'
c     in protein 'iprot', or 0 if it cannot be found
c
c     (C) Christopher Woods
c
c############################################################

      integer iprot,ires
      integer nres,i
      
      findResidue = 0
      
      if (iprot.le.0 .or. iprot.gt.MAXPROTEINS) return
      nres = NResidues(iprot)
      if (nres.le.0 .or. nres.gt.MAXRESIDUES) return
      
      do i=1,nres
        if (ResNum(iprot,i).eq.ires) then
          findResidue = i
          return
        endif
      enddo
      
      return
      end
