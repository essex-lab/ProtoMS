      
      integer function findSoluteAtom(isol,iatm,ires)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'solutes.inc'
      include 'enums.inc'
c#############################################################
c
c     This returns the zmatrix id of solute atom inam,
c     residue ires, in solute isol, or 0 if it doesn't
c     exist! -1, -2 or -3 are used to represent the 
c     special first three dummy atoms (residue name
c     'dum', atom names 'dm1', 'dm2', 'dm3')
c
c     (C) Christopher Woods
c
c#############################################################

      integer isol,i,nats
      character*4 iatm,ires
      
      findSoluteAtom = 0
      
      if (isol.le.0 .or. isol.gt.NSolutes) return
      nats = NSolAtoms(isol)
      if (nats.le.0) return

c     first return -1 -> -3 if it is a dummy...
      if (ires.eq.'dum') then
        if (iatm.eq."dm1") then
          findSoluteAtom = -1
          return
        else if (iatm.eq."dm2") then
          findSoluteAtom = -2
          return
        else if (iatm.eq."dm3") then
          findSoluteAtom = -3
          return
        endif
      endif
      
      do i=1,nats
        if (SolResNam(isol,SolResidue(isol,i)).eq.ires) then
          if (SolAtomNam(isol,i).eq.iatm) then
            findSoluteAtom = i
            return
          endif
        endif
      enddo

      write(printstring,*) "Could not find solute atom (",iatm,")-(",ires,") in solute ",SoluteName(isol)
      call printLine(WARNING,printstring)
      call printLine(WARNING,"Check your coordinate file to see that it has all the atoms that in the template")

      return
      end
