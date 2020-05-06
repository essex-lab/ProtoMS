
      subroutine processGeometryLine(check)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'solvents.inc'
      include 'splitter.inc'
      include 'constants.inc'
      include 'enums.inc'
c#########################################################
c
c     Process the 'geometry' input line - if check
c     is true then only print out what we would do.
c
c     The line has been split into words that are held
c     in the 'words' array
c
c     (C) Christopher Woods, January 2005
c
c#########################################################

      logical check
      integer ityp,imol,ires,istrt,igeom
      integer lineLength,ilgth
      integer getInt
      double precision getDouble
      character*4 inam,iresnam
      double precision ival
      
      if (nwords.lt.5) then
        call printLine(WARNING,"Could not interpret geometry line")
        call join(1,nwords,printstring)
        call printLine(WARNING,printstring)
        return
      endif
            
c     get the id of the molecule involved; format is
c     geometry solute isol ...
c     geometry solvent isvn ...
c     geometry protein iprot ires ...

      call lowercase(words(2))
      
      if (words(2).eq.'solute') then
        ityp = SOLUTE
        imol = getInt(words(3),0)
        if (imol.lt.1 .or. imol.gt.NSolutes) then
          write(printstring,*) "Requesting invalid solute! ",imol
          call printLine(WARNING,printstring)
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)
          return
        endif
      
        ilgth = lineLength(SoluteName(imol))
        write(printstring,*) "Modifying the geometry of solute ",imol," (",SoluteName(imol)(1:ilgth),")"
        call printLine(INFO,printstring)
        istrt = 4
      else if (words(2).eq.'protein') then
        ityp = PROTEIN
        imol = getInt(words(3),0)
        ires = getInt(words(4),0)
        
        if (imol.lt.1 .or. imol.gt.NProteins) then
          write(printstring,*) "Requesting invalid protein! ",imol
          call printLine(WARNING,printstring)
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)
          return
        endif

        if (ires.lt.1 .or. ires.gt.NResidues(imol)) then
          write(printstring,*) "Requesting invalid residue! ",ires
          call printLine(WARNING,printstring)
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)
          return
        endif
        
        ilgth = lineLength(ProteinName(imol))
        write(printstring,*) "Modifying the geometry of residue ",ResNam(imol,ires),"-",ResNum(imol,ires),
     .                        " in protein ",imol," (",ProteinName(imol)(1:ilgth),")"
        call printLine(INFO,printstring)
        istrt = 5
      else if (words(2).eq.'solvent') then
        ityp = SOLVENT
        imol = getInt(words(3),0)
        
        if (imol.lt.1 .or. imol.gt.NSolvents) then
          write(printstring,*) "Requesting invalid solvent! ",imol
          call printLine(WARNING,printstring)
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)
          return
        endif
        
        write(printstring,*) "Modifying the geometry of solvent ",imol," (",SvnName(imol),")"
        call printLine(INFO,printstring)
        istrt = 4
      else
        call printLine(WARNING,"Could not interpret the type of molecule whose geometry shold be changed")
        call join(1,nwords,printstring)
        call printLine(WARNING,printstring)
        return
      endif

c     now read the rest of the line; format is
c     geometry ... bond|angle|dihedral zmat            # resets bond/ang/dih to zmat values
c     geometry ... bond iatm 2.5                       # set the zmatrix bond of iatm to 2.5 A
c     geometry ... angle|dihedral iatm 30.0            # set the zmatrix ang/dih of iatm to 30 degrees
c     geometry ... bond|angle equilibrium              # resets bond/ang to equilibrium values
c                                                      # (taken from harmonic force constants)

c     make sure that the line is long enough
      if (nwords.lt.istrt+1) then
        call printLine(WARNING,"Geometry line is not long enough to interpret!")
        call join(1,nwords,printstring)
        call printLine(WARNING,printstring)
        return
      endif

      call lowercase(words(istrt))
      if (words(istrt).eq.'bond') then
        igeom = BONDS      
      else if (words(istrt).eq.'angle') then
        igeom = ANGLES      
      else if (words(istrt).eq.'dihedral') then
        igeom = DIHEDRALS
      else
        call printLine(WARNING,"Could not interpret what geometry should be changed!")
        call join(1,nwords,printstring)
        call printLine(WARNING,printstring)
        return
      endif

      call lowercase(words(istrt+1))
      
      if (words(istrt+1).eq.'zmat') then
        if (.not.check) call setZMatrixGeometry(ityp,imol,ires,igeom)
      else if (words(istrt+1).eq.'equilibrium') then
        if (igeom.eq.DIHEDRALS) then
          call printLine(WARNING,"Cannot set equilibrium geometry of dihedrals!")
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)
          return
        endif
        if (.not.check) call setEquilibriumGeometry(ityp,imol,ires,igeom)
      else
        istrt = istrt + 1
c       get the atom involved
        if (ityp.eq.PROTEIN) then
c         just need the atom name
          if (nwords.lt.istrt+1) then
            call printLine(WARNING,"Geometry line is not long enough to interpret!")
            call join(1,nwords,printstring)
            call printLine(WARNING,printstring)
            return
          endif
          
          call lowercase(words(istrt))
          inam = words(istrt)(1:4)
          
          ival = getDouble(words(istrt+1),ZERO)
          
          if (.not.check) call setProteinGeometry(imol,ires,inam,igeom,ival)
        else if (ityp.eq.SOLUTE) then
c         get the atom name and residue name
          if (nwords.lt.istrt+2) then
            call printLine(WARNING,"Geometry line is not long enough to interpret!")
            call join(1,nwords,printstring)
            call printLine(WARNING,printstring)
            return
          endif
          
          call lowercase(words(istrt))
          call lowercase(words(istrt+1))
          inam = words(istrt)(1:4)
          iresnam = words(istrt+1)(1:4)
          
          ival = getDouble(words(istrt+2),ZERO)
          
          if (.not.check) call setSoluteGeometry(imol,iresnam,inam,igeom,ival)
        else if (ityp.eq.SOLVENT) then
          call printLine(WARNING,"Cannot set the internal geometry of solvent molecules!")
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)
          return
        else
          call printLine(WARNING,"Cannot interpret type of molecule from geometry line")
          call join(1,nwords,printstring)
          call printLine(WARNING,printstring)
          return
        endif
      endif
              
      return
      end
