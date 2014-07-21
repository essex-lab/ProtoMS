      integer function findSolventAtom(isvn,inam)
      include 'dimensions.inc'
      include 'solvents.inc'
c###################################################
c
c     this returns the index of atom name inam in 
c     solvent isvn, or 0 if it doesn't exist
c
c     (C) Christopher Woods   5 March 2002
c
c###################################################

      integer isvn,i,nats
      character*4 inam
      
      findSolventAtom = 0
      if (isvn.le.0 .or. isvn.gt.NSolvents) return
      nats = NSvnAtoms(isvn)
      if (nats.le.0) return
      
      do i=1,nats
        if (SvnAtomName(isvn,i).eq.inam) then
          findSolventAtom = i
          return
        endif
      enddo

      return
      end
