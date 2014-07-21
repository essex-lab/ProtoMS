      
      integer function findAtom(iprot,ires,inam)
      include 'dimensions.inc'
      include 'proteins.inc'
c#######################################################
c
c     This function will attempt to find an atom with
c     atomname inam, in residue ires of iprotein...
c     It will return its zmatrix number if it finds it,
c     else 0 if it doesn't...
c
c     (C) Christopher Woods, October 2004
c
c#######################################################

      integer iprot,ires
      character*4 inam
      integer i,nat

      findAtom = 0

      if (iprot.le.0 .or. iprot.gt.NProteins) return
      if (ires.le.0 .or. ires.gt.NResidues(iprot)) return

c     first eliminate BBAtoms
      do i=1,4
        if (inam.eq.BBNam(iprot,ires,i)) then
          findAtom = BBZid(iprot,ires,i)
          return
        endif
      enddo
         
      nat = NSCAtoms(iprot,ires)
      if (nat.gt.0) then
        do i=1,nat
            if (SCNam(iprot,ires,i).eq.inam) then
              findAtom = SCZid(iprot,ires,i)
              return
            endif
        enddo
      endif
      
      return 
      
      end
