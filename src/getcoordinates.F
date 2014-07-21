
      subroutine getCoordinates(ityp,imol,iatm,ic,icb,icf)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'enums.inc'
      include 'solutes.inc'
      include 'proteins.inc'
      include 'printer.inc'
c################################################################################
c
c     This routine finds the coordinates of the atom specified in 
c     ityp,imol,ires,iatm (e.g. PROTEIN,1,102 = atom 102 of protein 1)
c     and returns them in 'ic', and the forwards and backwards in icf and icb.
c
c     If iatm=0, then the centre of the molecule is returned
c
c     (C) Christopher Woods  2 April 2003
c
c################################################################################

      integer ityp,imol,iatm,i
      double precision ic(3),icb(3),icf(3)
      
      if (ityp.eq.PROTEIN) then
c       we are looking for a protein atom!
        if (imol.lt.1 .or. imol.gt.NProteins) then
          write(printstring,*) "Cannot find protein ",imol
          call closeProgram(-1,printstring)
        endif
c       if iatm < 0, then return the centre of the protein residue numbered '-iatm'
        if (iatm.lt.0) then
          if (-iatm.gt.NResidues(imol)) then
            write(printstring,*) "Could not find residue ",-iatm," in protein ",imol
            call closeProgram(-1,printstring)
          endif
c         copy the coordinates of the centre of the residue to ic, icb and icf          
          do i=1,3
            ic(i) = ResCentre(imol,-iatm,i)
            icb(i) = ic(i)
            icf(i) = ic(i)
          enddo
        else if (iatm.gt.0) then
          if (iatm.gt.NZAtoms(imol)) then
            write(printstring,*) "Cannot find atom ",iatm," in protein ",imol
            call closeProgram(-1,printstring)
          endif
c         copy the coordinates of the protein atom to ic, icb and icf
          do i=1,3
            ic(i) = ZCoords(imol,iatm,i)
            icb(i) = ic(i)
            icf(i) = ic(i)
          enddo
        else
          write(printstring,*) "Cannot ask for atom number 0 in a protein!"
          call closeProgram(-1,printstring)
        endif      
      else if (ityp.eq.SOLUTE) then
c       we are now looking for a solute atom!
        if (imol.lt.1 .or. imol.gt.NSolutes) then
          write(printstring,*) "Cannot find solute ",imol
          call closeProgram(-1,printstring)
        endif
c       if iatm=0, then return the centre of the solute
        if (iatm.eq.0) then
          do i=1,3
            ic(i) = SolCentre(imol,i)
            icb(i) = ic(i)
            icf(i) = ic(i)
          enddo
        else if (iatm.gt.0) then
c         otherwise return the coordinates of the solute atom iatm
          if (iatm.gt.NSolAtoms(imol)) then
            write(printstring,*) "Cannot find atom ",iatm," in solute ",imol
            call closeProgram(-1,printstring)
          endif
c         if the solute has variable geometry, then return the forwards and backwards
c         coordinates correctly
          if (NSolVariable(imol).gt.0) then
            do i=1,3
              ic(i) = ZSolCoords(imol,iatm,i)
              icb(i) = ZSolCoords(NSolBackwards(imol),iatm,i)
              icf(i) = ZSolCoords(NSolForwards(imol),iatm,i)
            enddo
          else
c           otherwise just return the coordinates in the reference solute
            do i=1,3
              ic(i) = ZSolCoords(imol,iatm,i)
              icb(i) = ic(i)
              icf(i) = ic(i)
            enddo
          endif
        else
          write(printstring,*) "Cannot ask for a negative numbered solute atom!"
          call closeProgram(-1,printstring)
        endif
      else if (ityp.eq.SOLVENT) then
c       now we are trying to get a solvent atom!
        if (imol.lt.1 .or. imol.gt.NSolvents) then
          write(printstring,*) "Cannot find solvent molecule ",imol
          call closeProgram(-1,printstring)
        endif
c       if iatm=0 then return the centre of the solvent
        if (iatm.eq.0) then
          do i=1,3      
            ic(i) = SvnCentre(imol,i)
            icb(i) = ic(i)
            icf(i) = ic(i)
          enddo
        else if (iatm.gt.0) then
c         return the coordinates of solvent atom iatm in solvent molecule imol
          if (iatm.gt.NSvnAtoms(imol)) then
            write(printstring,*) "Could not find atom ",iatm," in solvent ",imol
            call closeProgram(-1,printstring)
          endif
          do i=1,3
            ic(i) = ZSvnCoords(imol,iatm,i)
            icb(i) = ic(i)
            icf(i) = ic(i)
          enddo
        else
          write(printstring,*) "Cannot ask for a negative numbered solvent atom!"
          call closeProgram(-1,printstring)
        endif
      else
        write(printstring,*) "Could not identify molecule type ",ityp
        call closeProgram(-1,printstring)
      endif
      
      return
      end
