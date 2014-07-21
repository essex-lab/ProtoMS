      integer function findChainAtom(itmpl,iatm)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'enums.inc'
c########################################################
c
c     This function returns the index of the atom
c     with name iatm in the chain template itmpl,
c     or 0 if the atom in not in the template.
c
c     (C) Christopher Woods 
c
c########################################################

      integer itmpl,i
      character*4 iatm
      
      findChainAtom = 0
      
      if (itmpl.le.0 .or. itmpl.gt.NChainTemplates) then
        call printLine(WARNING,"Trying to find atoms in a non-existant chain template!")
        call printLine(WARNING,"This is a bug!")
        return
      endif
      
c     see if this atom is one of the bbatoms
      do i=1,4
        if (ChainBBAtom(itmpl,i).eq.iatm) then
          findChainAtom = -i
          return
        endif
      enddo      
      
      if (NChainAtoms(itmpl).le.0) return
      
      do i=1,NChainAtoms(itmpl)
        if (ChainAtm(itmpl,i,1).eq.iatm) then
          findChainAtom = i
          return
        endif
      enddo
      
      return
      end
