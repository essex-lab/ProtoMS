      subroutine printPDBResidue(iprot,ires,showdummies,printstd)
      include 'dimensions.inc'
      include 'proteins.inc'
c#######################################################
c
c     This subroutine prints out residue ires, from protein
c     iprotein, to unit iunit, via the pdb format...
c
c     (C) Christopher Woods  21 January 2002
c
c#######################################################

      integer iprot,ires
      integer i,natoms
      logical showdummies,printstd,take
      
C     first print the backbone atoms
      do i=1,4
        call printPDBAtom(BBNum(iprot,ires,i),BBNam(iprot,ires,i),ResNam(iprot,ires),
     .                    ResNum(iprot,ires),BBCoord(iprot,ires,i,1),BBCoord(iprot,ires,i,2),
     .                    BBCoord(iprot,ires,i,3),ZPar(iprot,BBZid(iprot,ires,i)),showdummies,printstd)
      enddo
           
C     now print out the sidechain atoms...
      natoms = NSCAtoms(iprot,ires)
      if (natoms.gt.0) then
        do i=1,natoms
           take = .true.
           if (printstd.and.(SCNum(iprot,ires,i).eq.0)) take = .false.
           if (take) 
     .      call printPDBAtom(SCNum(iprot,ires,i),SCNam(iprot,ires,i),
     .      ResNam(iprot,ires),ResNum(iprot,ires),SCCoord(iprot,ires,i,1),
     .      SCCoord(iprot,ires,i,2),SCCoord(iprot,ires,i,3),
     .      ZPar(iprot,SCZid(iprot,ires,i)),showdummies,printstd)
        enddo
      endif     

      return
      end
