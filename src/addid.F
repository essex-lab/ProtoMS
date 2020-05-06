      subroutine addid()
      include 'dimensions.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'solvents.inc'
      include 'printer.inc'      
      include 'splitter.inc'
      include 'enums.inc'
      include 'labels.inc'
c#############################################
c     
c     This routine sets an assigned ID number to 
c     an atom 
c     The format of an clearid line is :
c
c     add 1 SOLUTE 1 C19 3T4
c     --> Assign id 1 to the atom C19 in the residue named 3T4 of Solute 1
c     add 4 PROTEIN 1 CA 34
c     --> Assign id 4 to the atom CA in the residue number 34 (PDB number)
c     of protein 1
c     add 6 SOLVENT 447 OH
c     --> Assign id 6 to the atom OH in the solvent number 447
c     
c     
c     Julien Michel, February 2005
c
c############################################
      integer atomid,ityp,imol,zid,pdbnumber,resnumber
      character*4 resname
      character*4 atname
      integer findAtom
      integer findSoluteAtom
      integer findSolventAtom
      integer findSolResidue
      integer findResidue
      integer getInt

      if (nwords .gt. 7) then
        write(printstring,*) 'Could not interpret setid line !'
        call printLine(WARNING,printstring)
        call join(1,nwords,printstring)
        call printLine(WARNING,printstring)        
      endif
      atomid = getInt(words(3),0)
      call lowercase(words(4))      
      if ( words(4).eq.'solute') then
        ityp = SOLUTE
        imol = getInt(words(5),0)
        if (imol.lt.1 .or. imol.gt.NSolutes) then
          write(printstring,*) "Requesting invalid solute! ",imol
          call printLine(WARNING,printstring)
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)
          return
        endif
        call lowercase(words(6))
        call lowercase(words(7))
        atname = words(6)(1:4)
        resname = words(7)(1:4)
c       Get the residue number, since we don't know the number, we set it to 0
        resnumber = findSolResidue(imol,0,resname)
        zid = findSoluteAtom(imol,atname,resname)
      else if ( words(4) .eq. 'protein') then
        ityp = PROTEIN
        imol = getInt(words(5),0)
        if (imol.lt.1 .or. imol.gt.NProteins) then
          write(printstring,*) "Requesting invalid protein! ",imol
          call printLine(WARNING,printstring)
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)
          return
        endif
        call lowercase(words(6))
        atname = words(6)(1:4)
c       Should be The PDB residue number
        pdbnumber = getInt(words(7),0)
        resnumber = findResidue(imol,pdbnumber)
        resname = ResNam(imol,resnumber)
        zid = findAtom(imol,resnumber,atname)
      else if ( words(4) .eq. 'solvent') then
        ityp = SOLVENT
        imol = getInt(words(5),0)
        
        if (imol.lt.1 .or. imol.gt.NSolvents) then
          write(printstring,*) "Requesting invalid solvent! ",imol
          call printLine(WARNING,printstring)
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)
          return
        endif
        call lowercase(words(6))
        atname = words(6)(1:4)
c       There is no number for a solvent residue, so assume 1
        resnumber = 1
        resname = SvnName(imol)
        zid = findSolventAtom(imol,atname)
      else
c        call printLine(WARNING,"Could not interpret the type of molecule for the id")
        call join(1,nwords,printstring)
        call printLine(WARNING,printstring)
        return  
      endif
      
      if (resnumber .eq.0) then
        write(printstring,6) resname,words(4),imol
        call printLine(WARNING,printstring)
        return
      endif
 6    Format("Could not find residue ",A4," in ",A8," number ",I4)
      
      if (zid .eq. 0) then
        write(printstring,7) atname,words(4),imol
        call printLine(WARNING,printstring)
        return
      endif
 7    Format("Could not find atom ",A4," in ",A8," number ",I4)
      
      if (LabelID(atomid) .ne. 0) then
        write(printstring,*) 'ID ',atomid,' is already used !'
        call printLine(WARNING,printstring)
        return
      endif

      write(printstring,10) atomid,atname,zid,words(4),imol,resnumber,resname
      call printLine(INFO,printstring)
 10   Format("Setting id ",I4," Atom ",A4," zid ",I4," ityp ",A8,
     &  " imol ",I4," resnumber ",I4," resname ",A4)

      LabelID(atomid) = zid
      LabelType(atomid) = ityp
      LabelName(atomid) = atname
      LabelMolecule(atomid) = imol
      LabelResidueName(atomid) = resname
      LabelResidueNumber(atomid) = resnumber
      
      return 
      end
