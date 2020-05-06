
      subroutine makeVariableBond(iprot,iat1,iat2,delta)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
c############################################################
c
C     this routine makes the bond between zmatrix atoms
c     at1 and at2 variable (assuming it *is* a zmatrix
c     defined bond!
c
c     (C) Christopher Woods   23 January 2002
c
c############################################################

      integer iprot,iat1,iat2
      double precision delta
      character*13 nameString,inam1,inam2
      integer ires,iatom

c     sort by higher number (since zmatricies are
c     defined in numerical order...)
      if (iat1.eq.iat2) return
      call sortZBond(iat1,iat2)
      
      ires = Zid2ResID(iprot,iat1,1)
      iatom = Zid2ResID(iprot,iat1,2)

      if (iatom.gt.0) then
c       we are using the sidechain
        if (ZMatSC(iprot,ires,iatom,1).eq.iat2) then
          ZMovSC(iprot,ires,iatom,1) = .true.
          ZDeltaSC(iprot,ires,iatom,1) = delta
        else
          inam1 = nameString(iprot,iat1)
          inam2 = nameString(iprot,iat2)
          write(printstring,*) "Trying to make a non-zmatrix-defined internal coordinate flexible...",
     .              inam1,inam2
          call printLine(WARNING,printstring)
        endif
      else
          inam1 = nameString(iprot,iat1)
          inam2 = nameString(iprot,iat2)
          write(printstring,*) "Trying to make a non-zmatrix-defined internal coordinate flexible...",
     .              inam1,inam2
          call printLine(WARNING,printstring)
      endif
      
      return
      
      end
