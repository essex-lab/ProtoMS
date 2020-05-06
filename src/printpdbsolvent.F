      subroutine printPDBSolvent(isvn,showdummies,printstd)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'printer.inc'
c############################################################
c
c     This routine prints the solvent as a pdb
c
c     (C) Christopher Woods   31 January 2002
c
c############################################################

      integer i,nats,isvn
      logical showdummies,printstd
      
      if (isvn.le.0 .or. isvn.gt.NSolvents) return
      nats = NSvnAtoms(isvn)
      if (nats.le.0) return
      
      do i=1,nats
        call printPDBAtom(i,SvnAtomName(isvn,i),SvnName(isvn),SvnRes(isvn),ZSvnCoords(isvn,i,1),
     .                    ZSvnCoords(isvn,i,2),ZSvnCoords(isvn,i,3),ZSvnPar(isvn,i),showdummies,printstd)
      enddo

c     also print the axes and centre
c      call printPDBAtom(i,"CTR",SvnName(isvn),isvn,SvnCentre(isvn,1),SvnCentre(isvn,2),SvnCentre(isvn,3),0)
c      call printPDBAtom(i,"MAJ",SvnName(isvn),isvn,SvnCentre(isvn,1)+SvnMajAxis(isvn,1),
c     .                                             SvnCentre(isvn,2)+SvnMajAxis(isvn,2),
c     .                                             SvnCentre(isvn,3)+SvnMajAxis(isvn,3),0)
c      call printPDBAtom(i,"MIN",SvnName(isvn),isvn,SvnCentre(isvn,1)+SvnMinAxis(isvn,1),
c     .                                             SvnCentre(isvn,2)+SvnMinAxis(isvn,2),
c     .                                             SvnCentre(isvn,3)+SvnMinAxis(isvn,3),0)
c      call printPDBAtom(i,"PER",SvnName(isvn),isvn,SvnCentre(isvn,1)+SvnPerAxis(isvn,1),
c     .                                             SvnCentre(isvn,2)+SvnPerAxis(isvn,2),
c     .                                             SvnCentre(isvn,3)+SvnPerAxis(isvn,3),0)
c      

c     print a 'TER' line
      call printLine(PDB,"TER")
      
      return
      end
