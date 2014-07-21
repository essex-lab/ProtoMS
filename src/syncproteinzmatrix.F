      subroutine syncProteinZmatrix(iprot)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'proteins.inc'
c###########################################################
c
c     This subroutine rebuilds the zmatrix of protein iprot
c     from the coordinates currently in ZCoords
c
c     (C) Christopher Woods, October 2004
c
c###########################################################

      integer iprot,nres,i,j,k,natoms
      double precision ibnd,iang,idih,calcBond,calcAngle,calcDihedral
      integer ilgth,lineLength
      
      nres = NResidues(iprot)
      if (nres.le.0 .or. nres.gt.MAXRESIDUES) then
        write(printstring,*) "Trying to build a protein with an insane number of residues!"
        call printLine(WARNING,printstring)
        return
      endif
      
      do i=1,nres
c       copy the backbone coordinates from the zmatrix array to the residue-centric array
        do j=1,4
          do k=1,3
            if (BBZid(iprot,i,j).eq.0) then
              ilgth = lineLength(ProteinName(iprot))
              write(printstring,*) "Serious problem with bbatom ",BBNam(iprot,i,j)," in ",
     .                  ProteinName(iprot)(1:ilgth)," ",ResNam(iprot,i),"-",ResNum(iprot,i)
              call printLine(WARNING,printstring)
              call printLine(WARNING,"This is a bug, and your simulation is now broken!")
            else
              BBCoord(iprot,i,j,k) = ZCoords(iprot,BBZid(iprot,i,j),k)                  
            endif
          enddo
        enddo

        natoms = NSCAtoms(iprot,i)
        do j=1,natoms
c         update the residue centric coordinates from the zmatrix coordinates
          do k=1,3
            if (SCZid(iprot,i,j).eq.0) then
              ilgth = lineLength(ProteinName(iprot))
              write(printstring,*) "Serious problem with atom ",SCNam(iprot,i,j)," in ",
     .                  ProteinName(iprot)(1:ilgth)," ",ResNam(iprot,i),"-",ResNum(iprot,i)
              call printLine(WARNING,printstring)
              call printLine(WARNING,"This is a bug, and your simulation is now broken!")
            else
              SCCoord(iprot,i,j,k) = ZCoords(iprot,SCZid(iprot,i,j),k)                  
            endif
          enddo
c         calculate the internal coordinates
          ibnd = calcBond(PROTEIN,iprot,SCZid(iprot,i,j),ZMatSC(iprot,i,j,1))
          iang = calcAngle(PROTEIN,iprot,SCZid(iprot,i,j),ZMatSC(iprot,i,j,1),ZMatSC(iprot,i,j,2))
          idih = calcDihedral(PROTEIN,iprot,SCZid(iprot,i,j),ZMatSC(iprot,i,j,1),
     .              ZMatSC(iprot,i,j,2),ZMatSC(iprot,i,j,3))
     
c         now assign these properly
          ZValSC(iprot,i,j,1) = ibnd
          ZValSC(iprot,i,j,2) = iang
          ZValSC(iprot,i,j,3) = idih
        enddo
      enddo
     
c     recalculate the centre and radius of each residue
      call getResidueCentreRadii(iprot)
c     now recalculate all of the bond lengths and angle/dihedral sizes
      do i=1,nres
        call updateResidueInternals(iprot,i)
      enddo
      return
     
      end
