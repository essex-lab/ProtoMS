      subroutine updateBond(iprot,ires,iresbnd)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'enums.inc'
c###########################################################################
c
c     This routine updates the bond parameters for bond iresbnd in 
c     residue ires, protein iprot. (only for bonds which are in the 
c     template itmpl)
c
c     (C) Christopher Woods  4 February 2002
c
c###########################################################################

      integer iprot,ires,iresbnd
      integer ibnd,zat1,zat2
      integer getBond
      
      if (iresbnd.le.0 .or. iresbnd.gt.NResBonds(iprot,ires)) return
      
c     get the bond number in the protein of this bond...
      ibnd = ResidueBond(iprot,ires,iresbnd)
      if (Bond(iprot,ibnd,3).lt.0) then
        write(printstring,*) "Attempting to convert a residue with variable ",
     .                        "bond parameters! - ",iprot,ires
        call closeProgram(-1,printstring)
        return
      endif
c     do not convert bonds that are not part of the template or backbone
c     (these are null energy bonds)
      if (BondDummy(iprot,ibnd)) return
      zat1 = Bond(iprot,ibnd,1)
      zat2 = Bond(iprot,ibnd,2)

      Bond(iprot,ibnd,3) = getBond(ZPar(iprot,zat1),ZPar(iprot,zat2))
      
      return
      
      end
