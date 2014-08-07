      subroutine writeRestart
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'simparams.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'stacks.inc'
c########################################################
c
c     This subroutine writes out the information
c     needed to perform a restart from this point
c
c     (C) Christopher Woods   25 November 2002
c
c########################################################

      integer i,j,k,ires,iat
      integer isvn,igcs
      character*4 inam
      character*4 iresnam
      double precision savenrg,inrgb,inrgf,ienew,totalEnergy
c     The restart file is just a list of all proteins,
c     then solutes, then solvent coordinates. The aim
c     still is that the next simulation loads up the same
c     input files, then the restart file is used to modify
c     the coordinates
c      ienew = totalEnergy(savenrg,inrgb,inrgf)
c        PRINT*, "ENERGY",ienew
c     don't do anything if nothing is loaded
      if (NProteins.le.0 .and. NSolutes.le.0 .and. NSolvents.le.0)
     .             return
c      PRINT*, "GCSCavity-nwrite",GCSCavity
c      write(printstring,10) "Simulation B factor" ,B
c      call printLine(DETAIL,printstring)
c     now write out the number of proteins,solutes,solvents
      write(printstring,10) "NPROTEINS= ",NProteins
      call printLine(RESTART,printstring)
      write(printstring,10) "NSOLUTES= ",NSolutes
      call printLine(RESTART,printstring)
      write(printstring,10) "NSOLVENTS= ",NSolvents
      call printLine(RESTART,printstring)
      write(printstring,10) "GCSOLUTES= ",NGCSolutes
      call printLine(RESTART,printstring)
 
10    format(a11,i8)

       write(printstring,91) GCSCavity
       call printLine(RESTART,printstring)
91     format('GCSCavity= ',i8)

c     if necessary, write out information about the solvent box
      if (boundary.eq.PERIODICBOUNDARIES) then
        write(printstring,72) BoxDimension(1),BoxDimension(2),BoxDimension(3)
        call printLine(RESTART,printstring)
        write(printstring,72) BoxOrigin(1),BoxOrigin(2),BoxOrigin(3)
        call printLine(RESTART,printstring)
72      format(f25.10,1x,f25.10,1x,f25.10)
      endif
      
c     now write out the coordinates of each protein
      do i=1,NProteins
        write(printstring,*) "Writing restart info for protein ",i,"..."
        call printLine(DETAIL,printstring)
        write(printstring,8) i,ProteinName(i)
8       format('Protein ',i8,' = ',a)
        call printLine(RESTART,printstring)
        write(printstring,12) NZAtoms(i)
        call printLine(RESTART,printstring)
12      format('NAtoms = ',i8)        
        
        do j=1,NZAtoms(i)
          ires = Zid2ResID(i,j,1)
          iat = Zid2ResID(i,j,2)
          if (iat.gt.0) then
            inam = SCNam(i,ires,iat)
          else if (iat.lt.0) then
            inam = BBNam(i,ires,-iat)
          else
            inam = "????"
          endif 
          if (ires.gt.0) then
            iresnam = ResNam(i,ires)
          else
            iresnam = "???"
          endif
          
          write(printstring,11) i,j,inam,iresnam,ires,
     .                                     (ZCoords(i,j,k),k=1,3)
11        format(i8,1x,i8,1x,a4,1x,a3,1x,i8,1x,3(f14.8,1x))
          call printLine(RESTART,printstring)
        enddo
      enddo

c     now write the coordinates of all of the solutes
      do i=1,NSolutes
        write(printstring,*) "Writing restart info for solute ",i,"..."
        call printLine(DETAIL,printstring)
        write(printstring,20) i,SoluteName(i)
20      format('Solute ',i8,' = ',a)
        call printLine(RESTART,printstring)
        write(printstring,21) NSolAtoms(i)
        call printLine(RESTART,printstring)
21      format('NAtoms = ',i8)
c       SG
        write(printstring,77) SoluteTheta(i)
        call printLine(RESTART,printstring)        
77      format('Theta = ',f8.4)
c       first write out the dummies / axes of the solute
        write(printstring,23) 'DUMMY',(SolDummy(i,1,j),j=1,3)
        call printLine(RESTART,printstring)
        write(printstring,23) 'DUMMY',(SolDummy(i,2,j),j=1,3)
        call printLine(RESTART,printstring)
        write(printstring,23) 'DUMMY',(SolDummy(i,3,j),j=1,3)
        call printLine(RESTART,printstring)
        write(printstring,23) 'GCENT',(SolGeomCentre(i,j),j=1,3)
        call printLine(RESTART,printstring)
c        write(printstring,23) 'MAJOR',(SolAxis(i,1,j),j=1,3)
c        call printLine(RESTART,printstring)
c        write(printstring,23) 'MINOR',(SolAxis(i,2,j),j=1,3)
c        call printLine(RESTART,printstring)
c        write(printstring,23) 'PERPD',(SolAxis(i,3,j),j=1,3)
c        call printLine(RESTART,printstring)
23      format(a5,1x,3(f14.8,1x))

        do j=1,NSolAtoms(i)
          write(printstring,22) i,j,SolAtomNam(i,j),
     .        SolResNam(i,SolResidue(i,j)),SolResidue(i,j),
     .        (ZSolCoords(i,j,k),k=1,3)
22        format(i8,1x,i8,1x,a4,1x,a3,1x,i8,1x,3(f14.8,1x))
          call printLine(RESTART,printstring)
        enddo
      enddo

c     now write the coordinates of all of the solvents
      if (NSolvents.gt.0) then
        write(printstring,*) "Writing restart info for the solvent..."
        call printLine(DETAIL,printstring)
        do i=1,NSolvents
          isvn = i

          write(printstring,30) isvn,SvnName(i)
30        format('Solvent ',i8,' = ',a4)
          call printLine(RESTART,printstring)
          write(printstring,31) NSvnAtoms(i)
31        format('NAtoms = ',i8)
          call printLine(RESTART,printstring)

c         now print out the axis information for the solvent
          write(printstring,33) 'CNTRE',(SvnCentre(i,j),j=1,3)
          call printLine(RESTART,printstring)
c          write(printstring,33) 'MAJOR',(SvnMajAxis(i,j),j=1,3)
c          call printLine(RESTART,printstring)
c          write(printstring,33) 'MINOR',(SvnMinAxis(i,j),j=1,3)
c          call printLine(RESTART,printstring)
c          write(printstring,33) 'PERPD',(SvnPerAxis(i,j),j=1,3)
c          call printLine(RESTART,printstring)
33        format(a5,1x,3(f14.8,1x))

          do j=1,NSvnAtoms(i)
            write(printstring,32) isvn,j,SvnAtomName(i,j),SvnName(i),
     .           (ZSvnCoords(i,j,k),k=1,3)
32          format(i8,1x,i8,1x,a4,1x,a3,1x,3(f14.8,1x))
            call printLine(RESTART,printstring)
          enddo
        enddo
      endif
      
c       now write the coordinates of all of the GC-solutes
      if (NGCSolutes.gt.0) then
        write(printstring,*) "Writing restart info for GC solutes..."
        call printLine(DETAIL,printstring)
        do i=1,NGCSolutes
          igcs = i

          write(printstring,40) igcs,GCSName(i)
40        format('GC-solute ',i8,' = ',a4)
          call printLine(RESTART,printstring)

          write(printstring,99) thetaparam(i)
99        format('Theta = ',g20.8)
          call printLine(RESTART,printstring)

          write(printstring,41) NGCSAtoms(i)
41        format('NAtoms = ',i8)
          call printLine(RESTART,printstring)

c         now print out the axis information for the solvent
          write(printstring,43) 'CNTRE',(GCSCentre(i,j),j=1,3)
          call printLine(RESTART,printstring)
43        format(a5,1x,3(f14.8,1x))

          do j=1,NGCSAtoms(i)
            write(printstring,42) igcs,j,GCSAtomName(i,j),GCSName(i),
     .           (ZGCSCoords(i,j,k),k=1,3)
42          format(i8,1x,i8,1x,a4,1x,a3,1x,3(f14.8,1x))
            call printLine(RESTART,printstring)
          enddo
        enddo
      endif
      

      return
      end
