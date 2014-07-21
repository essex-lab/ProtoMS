      subroutine sumCharge(iprot)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'constants.inc'
      include 'proteins.inc'
c############################################################
c
c     This subroutine sums the charge on the whole protein
c     and prints it out...
c
c     (C) Christopher Woods   30 January 2002
c
c############################################################

      integer iprot,i,j,l,nres
      double precision ichg(MAXRESIDUES),totalchg
      double precision sumResidueCharge

      if (iprot.le.0 .or. iprot.gt.NProteins) return
      
      nres = NResidues(iprot)
      
      if (nres.le.0) return
      
      totalchg = ZERO
      do i=1,nres
        ichg(i) = sumResidueCharge(iprot,i)
        write(printstring,*) 'Residue ',i,' charge ',ichg(i)
        call printLine(DETAIL,printstring)
        totalchg = totalchg + ichg(i)
      enddo

      write(printstring,10) iprot,totalchg
      call prettyLine(INFO,printstring)
10    format("Protein ",i2," has a total charge of ",f12.5)
      
      do i=1,nres,8
        j = i+8
        if (j.gt.nres) j = nres
        write(printstring,11) (ResNam(iprot,l),ResNum(iprot,l),ichg(l),l=i,j)
        call printLine(DETAIL,printstring)
      enddo
11    format(30(" ",a,"(",i3,")[",f6.3,"]"))

      TotalCharge(iprot) = totalchg

      return
      end
