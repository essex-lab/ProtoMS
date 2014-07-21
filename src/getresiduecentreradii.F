      subroutine getResidueCentreRadii(iprot)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'proteins.inc'
c##############################################################
c
c     This routine calculates all of the residue centres and
c     radii for protein iprot
c
c     (C) Christopher Woods   6 March 2002
c
c##############################################################

      integer iprot,i,nres
      
      if (iprot.gt.NProteins .or. iprot.le.0) return
      
      nres = NResidues(iprot)
      if (nres.le.0 .or. nres.gt.MAXRESIDUES) return
      
c      write(printstring,*) "Getting the centres and radii of residues in protein ",iprot
c      call printLine(INFO,printstring)
      
      do i=1,nres
        call getResidueCentreRadius(iprot,i)
      enddo
      return
      end
