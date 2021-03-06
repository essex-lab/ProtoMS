      subroutine addAcid(iprot,ires)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
c##############################################################
c
c     This adds residue iprot,ires to the acid array...
c
c     (C) Christopher Woods   4 February 2002
c
c##############################################################

      integer iprot,ires
      
      NAcidicRes(iprot) = NAcidicRes(iprot) + 1
      if (NAcidicRes(iprot).gt.MAXACIDICRES) then
        write(printstring,*) "Exceeded maximum number of acidic residues, ",MAXACIDICRES
        call printLine(FATAL,printstring)
        write(printstring,*)"Increase MAXACIDICRES and recompile"
        call closeProgram(-1,printstring)
        NAcidicRes(iprot) = NAcidicRes(iprot) - 1
        return
      endif
      
      AcidicRes(iprot,NAcidicRes(iprot)) = ires
      
      return
      end
