
      subroutine makeVariableAngle(iprot,iat1,iat2,iat3,delta)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
c############################################################
c
C     this routine makes the angle between zmatrix atoms
c     at1, at2 and at3 variable (assuming it *is* a zmatrix
c     defined angle!
c
c     (C) Christopher Woods   23 January 2002
c
c############################################################

      integer iprot,iat1,iat2,iat3
      integer ires,iatom
      character*13 nameString,inam1,inam2,inam3
      double precision delta

c     sort by higher number (since zmatricies are
c     defined in numerical order...)
      if (iat1.eq.iat3) return
      call sortZAngle(iat1,iat2,iat3)
      
      ires = Zid2ResID(iprot,iat1,1)
      iatom = Zid2ResID(iprot,iat1,2)

      if (iatom.gt.0) then
c       we are using the sidechain
        if (ZMatSC(iprot,ires,iatom,1).eq.iat2 .and.
     .        ZMatSC(iprot,ires,iatom,2).eq.iat3) then
            ZMovSC(iprot,ires,iatom,2) = .true.
            ZDeltaSC(iprot,ires,iatom,2) = delta
        else
          inam1 = nameString(iprot,iat1)
          inam2 = nameString(iprot,iat2)
          inam3 = nameString(iprot,iat3)
          write(printstring,*) "Trying to make a non-zmatrix-defined internal coordinate flexible...",
     .              inam1,inam2,inam3
          call printLine(WARNING,printstring)
        endif
      else
c       its a backbone, so cannot be defined by the zmatrix
          inam1 = nameString(iprot,iat1)
          inam2 = nameString(iprot,iat2)
          inam3 = nameString(iprot,iat3)
          write(printstring,*) "Trying to make a non-zmatrix-defined internal coordinate flexible...",
     .              inam1,inam2,inam3
        call printLine(WARNING,printstring)
      endif
      
      return
      
      end
