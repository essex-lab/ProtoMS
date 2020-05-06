      subroutine updateUreyBradley(iprot,ires,iresuby)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'enums.inc'
c###########################################################################
c
c     This routine updates the Urey-Bradley parameters for uby iresuby in 
c     residue ires, protein iprot.
c
c     (C) Christopher Woods
c
c###########################################################################

      integer iprot,ires,iresuby
      integer iuby,zat1,zat2,zat3
      integer findUreyBradleyParam
      
      if (iresuby.le.0 .or. iresuby.gt.NResUreyBradleys(iprot,ires)) return
      
c     get the uby number in the protein of this bond...
      iuby = ResidueUreyBradley(iprot,ires,iresuby)
      if (UreyBradley(iprot,iuby,3).lt.0) then
        write(printstring,*) "Attempting to convert a residue with variable ",
     .                        "Urey-Bradley parameters! - ",iprot,ires
        call closeProgram(-1,printstring)
        return
      endif

      zat1 = UreyBradley(iprot,iuby,1)
      zat2 = UreyBradley(iprot,iuby,2)
      zat3 = UreyBradley(iprot,iuby,3)

      UreyBradley(iprot,iuby,3) = findUreyBradleyParam(ZPar(iprot,zat1),ZPar(iprot,zat2),ZPar(iprot,zat3))
      
      return
      
      end
