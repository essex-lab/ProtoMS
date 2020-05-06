      subroutine printPDBProtein(iprotein,showdummies,printstd)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'proteins.inc'
c#######################################################
c
c     This subroutine prints the pdb version of protein
c     iprotein to either stdout, or filename
c
c     (C) Christopher Woods 21 January 2002
c
c#######################################################

      integer iprotein   
      integer nres,i
      logical showdummies,printstd

      if ((iprotein.gt.MAXPROTEINS).or.(iprotein.lt.1)) then
        write(printstring,*) "Requesting the printing of an invalid protein! ",iprotein
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is probably a programming bug!")
        return
      endif
      
      if (NResidues(iprotein).le.0) then
        write(printstring,*) "Trying to print a protein that contains no residues!"
        call printLine(WARNING,printstring)
        return
      endif

c     print all residues...
      nres = NResidues(iprotein)
      do i=1,nres
        call printPDBResidue(iprotein,i,showdummies,printstd)
      enddo

c     print a 'TER' line
      call printLine(PDB,"TER")
      
      end
