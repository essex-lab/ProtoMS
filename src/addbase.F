      subroutine addBase(iprot,ires)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
c##############################################################
c
c     This adds residue iprot,ires to the base array
c
c     (C) Christopher Woods  4 February 2002
c
c##############################################################

      integer iprot,ires
      
      NBasicRes(iprot) = NBasicRes(iprot) + 1
      if (NBasicRes(iprot).gt.MAXBASICRES) then
        write(printstring,*) "Exceeded maximum number of basic residues, ",MAXBASICRES
        call printLine(FATAL,printstring)
        write(printstring,*)"Increase MAXBASICRES and recompile"
        call closeProgram(-1,printstring)
        NBasicRes(iprot) = NBasicRes(iprot) - 1
        return
      endif
      
      BasicRes(iprot,NBasicRes(iprot)) = ires
      
      return
      end
