      integer function findMolTemplateAtom(itmp,inam,ires)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
c#########################################################
c
c     This function returns the template number of atm
c     inam, residue ires in molecule template itmp, or 0 if it doesn't exist
c
c     (C) Christopher Woods
c
c##########################################################

      integer itmp,i,nats
      character*4 inam
      character*4 ires
      
      findMolTemplateAtom = 0
      if (itmp.le.0 .or. itmp.gt.NMolTemplates) then
        write(printstring,*) "Trying to access impossible template (",itmp,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is probably a sign of a programming bug!")
        return
      endif
      
      nats = NMolTemplateAtms(itmp)
      if (nats.le.0) return
      
      do i=1,nats
        if (MolTemplateAtmName(itmp,i).eq.inam
     .             .and. MolTemplateAtmRes(itmp,i).eq.ires) then
          findMolTemplateAtom = i
          return
        endif
      enddo

      return
      
      end
