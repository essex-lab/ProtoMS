      
      subroutine readSolTemplateAtomLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'splitter.inc'
      include 'enums.inc'
c###################################################################
c
c     This routine reads in the atoms from a molecule template line
c
c     (C) Christopher Woods 
c
c###################################################################

      character*(*) line
      integer itmpl
      character*4 iatm1,ires1,bndatm,bndres,angatm,angres,dihatm,dihres
      integer findMolTemplateAtom
      integer iclj0,iclj1,iatm,istrt
      integer readInt
      
c     format is 'atom iatm ires iclj0 iclj1 bndatm bndres angatm angres dihatm dihres'
      if (nwords.ge.11) then
        iclj0 = readInt(words(4))
        iclj1 = readInt(words(5))
        istrt = 6
      else if (nwords.eq.10) then
c       format is 'atom iatm ires iclj bndatm bndres angatm angres dihatm dihres', ipar0=ipar1
        iclj0 = readInt(words(4))
        iclj1 = iclj0
        istrt = 5
      else if (nwords.eq.9) then
c       format is 'atom iatm ires bndatm bndres angatm angres dihatm dihres', ipar0=ipar1=0
        iclj0 = 0
        iclj1 = 0      
        istrt = 4
      else
        write(printstring,*) "Incorrectly formatted solute atom template line?"
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
        return
      endif
      
c     read in the line
      iatm1 = words(2)(1:4)
      ires1 = words(3)(1:4)
      bndatm = words(istrt)(1:4)
      bndres = words(istrt+1)(1:4)
      angatm = words(istrt+2)(1:4)
      angres = words(istrt+3)(1:4)
      dihatm = words(istrt+4)(1:4)
      dihres = words(istrt+5)(1:4)        
      
c     make sure that values of less than 0 are interpreted as dummies (0)
      if (iclj0.lt.0 .or. iclj0.gt.MAXCLJ) iclj0 = 0
      if (iclj1.lt.0 .or. iclj1.gt.MAXCLJ) iclj1 = 0

c     see if the atom already exists
      iatm = findMolTemplateAtom(itmpl,iatm1,ires1)
      if (iatm.ne.0) then
        write(printstring,*) "Overwriting information for solute atom (",iatm1,",",ires1,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
      else
        iatm = NMolTemplateAtms(itmpl) + 1
      
        if (iatm.gt.MAXMOLTEMPLATEATMS) then
          write(printstring,*) "Exceeded maximum number of solute template atoms (",
     .                         MAXMOLTEMPLATEATMS,")"
          call printLine(WARNING,printstring)
          call printLine(WARNING,"Increase MAXMOLTEMPLATEATMS and recompile")
          return
        endif
        NMolTemplateAtms(itmpl) = iatm
      endif
      
      MolTemplateAtmName(itmpl,iatm) = iatm1
      MolTemplateAtmRes(itmpl,iatm) = ires1
      MolTemplateAtmCLJ(itmpl,iatm,1) = iclj0
      MolTemplateAtmCLJ(itmpl,iatm,2) = iclj1
      MolTemplateZNam(itmpl,iatm,1) = bndatm
      MolTemplateZRes(itmpl,iatm,1) = bndres
      MolTemplateZNam(itmpl,iatm,2) = angatm
      MolTemplateZRes(itmpl,iatm,2) = angres
      MolTemplateZNam(itmpl,iatm,3) = dihatm
      MolTemplateZRes(itmpl,iatm,3) = dihres
      
      return
      end
