      subroutine printPDBGCSolute(igcs,showdummies,printstd)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'printer.inc'
      include 'stacks.inc'
c############################################################
c
c     This routine prints the solvent as a pdb
c
c     (C) Christopher Woods   31 January 2002
c
c############################################################

      integer i,nats,igcs
      logical showdummies,getOption,printstd
      
      if (igcs.le.0 .or. igcs.gt.NGCSolutes) return
      if (getOption('jaws2')) then
      goto 9
      endif
      if (thetaparam(igcs).lt.0.950) return
9     nats = NGCSAtoms(igcs)
      if (nats.le.0) return
      
      do i=1,nats
        call printPDBAtom(i,GCSAtomName(igcs,i),GCSName(igcs),GCSRes(igcs),ZGCSCoords(igcs,i,1),
     .                    ZGCSCoords(igcs,i,2),ZGCSCoords(igcs,i,3),ZGCSPar(igcs,i),showdummies,printstd)
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
