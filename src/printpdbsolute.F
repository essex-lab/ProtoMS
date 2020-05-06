      subroutine printPDBSolute(isol,showdummies,printstd)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'solutes.inc'
c#####################################################################
c
c     This subroutine prints out the solute as a pdb
c
c     (C) Christopher Woods   31 January 2002
c
c#####################################################################
      
      integer i,nats,isol
      logical showdummies,printstd
      
      if (isol.le.0 .or. isol.gt.NSolutes) return
      nats = NSolAtoms(isol)
      if (nats.le.0) return
 
c     print the first three dummy atoms if 'showdummies' is true
      if (showdummies) then
        call printPDBAtom(1,"dm1","dum",1,SolDummy(isol,1,1),SolDummy(isol,1,2),SolDummy(isol,1,3),0,.true.)
        call printPDBAtom(1,"dm2","dum",1,SolDummy(isol,2,1),SolDummy(isol,2,2),SolDummy(isol,2,3),0,.true.)
        call printPDBAtom(1,"dm3","dum",1,SolDummy(isol,3,1),SolDummy(isol,3,2),SolDummy(isol,3,3),0,.true.)
      endif
                   
      do i=1,nats
        call printPDBAtom(SolAtomNum(isol,i),SolAtomNam(isol,i),
     .                    SolResNam(isol,SolResidue(isol,i)),SolResNum(isol,SolResidue(isol,i)),
     .                    ZSolCoords(isol,i,1),ZSolCoords(isol,i,2), 
     .                    ZSolCoords(isol,i,3),ZSolPar(isol,i),showdummies,printstd)
      enddo

c     print a 'TER' line
      call printLine(PDB,"TER")
      
      return
      end
