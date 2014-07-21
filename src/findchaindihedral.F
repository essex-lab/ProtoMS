      integer function findChainDihedral(at1,at2,at3,at4,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c#######################################################
c
c      This function returns the index of the dihedral
c      between atoms at1-at2-at3-at4 in the chain template
c      itmpl, or 0 if there is no such dihedral
c
c      (C) Christopher Woods
c
c#######################################################

      integer itmpl,i
      character*4 at1,at2,at3,at4
      
      findChainDihedral = 0
      
      if (NChainDihedrals(itmpl).le.0) return
      
      call sortTemplateDihedral(at1,at2,at3,at4)
      
      do i=1,NChainDihedrals(itmpl)
        if (at1.eq.ChainDihName(itmpl,i,1) .and.
     .      at2.eq.ChainDihName(itmpl,i,2) .and.
     .      at3.eq.ChainDihName(itmpl,i,3) .and.
     .      at4.eq.ChainDihName(itmpl,i,4)) then
          findChainDihedral = i
          return
        endif
      enddo
      
      return
      end
