      subroutine addProteinUreyBradley(iprot,iat1,iat2,iat3,dummy,ipar0,ipar1)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'proteins.inc'
c######################################################
c
c     This function adds a urey-bradley between atoms with zmatrix
c     id's of zat1 and zat2...
c
c     (C) Christopher Woods
c
c######################################################

      integer iprot,zat1,zat2,zat3,ires1,ires2,ires3,nubys,iparam
      integer iat1,iat2,iat3,ipar0,ipar1
      logical dummy
      integer getUreyBradley,setDirectUreyBradley
      logical hasUreyBradley
      
      zat1 = iat1
      zat2 = iat2
      zat3 = iat3
      
c     sort the atoms
      call sortZUreyBradley(zat1,zat2,zat3)
      
      ires1 = Zid2ResID(iprot,zat1,1)
      ires2 = Zid2ResID(iprot,zat2,1)
      ires3 = Zid2ResID(iprot,zat3,1)

c     make sure that the uby hasn't already been added...
      if (hasUreyBradley(iprot,zat1,zat2,zat3)) then
        return
      endif
      
c     now add the uby...
      NUreyBradleys(iprot) = NUreyBradleys(iprot) + 1
      nubys = NUreyBradleys(iprot)
      if (nubys.gt.MAXUREYBRADLEYS) then
        write(printstring,*) "Exceeded maximum number of protein Urey-Bradley terms! (",
     .                        MAXUREYBRADLEYS,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Increase MAXUREYBRADLEYS and recompile")
        NUreyBradleys(iprot) = NUreyBradleys(iprot) - 1
        return
      endif
      UreyBradley(iprot,nubys,1) = zat1
      UreyBradley(iprot,nubys,2) = zat2
      UreyBradley(iprot,nubys,3) = zat3

      if (dummy) then
c       this is a dummy Urey Bradley (what's the point of that?)
        iparam = 0
      else if (ipar0.ge.0 .and. ipar1.ge.0) then
        iparam = setDirectUreyBradley(ipar0,ipar1)
      else
        iparam = getUreyBradley(ZPar(iprot,zat1),ZPar(iprot,zat2),ZPar(iprot,zat3)) 
      endif
      UreyBradley(iprot,NUreyBradleys(iprot),4) = iparam

c     now add the bond to the residue
      NResUreyBradleys(iprot,ires1) = NResUreyBradleys(iprot,ires1) + 1
      if (NResUreyBradleys(iprot,ires1).gt.MAXUREYBRADLEYSPERRESIDUE) then
        write(printstring,*) "Exceeded maximum number of Urey-Bradley terms in residue ",ires1
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Increase MAXUREYBRADLEYSPERRESIDUE and recompile")
        NResUreyBradleys(iprot,ires1) = NResUreyBradleys(iprot,ires1) - 1
        return
      endif
      ResidueUreyBradley(iprot,ires1,NResUreyBradleys(iprot,ires1)) = NUreyBradleys(iprot)
      if (ires1.ne.ires2) then
        NResUreyBradleys(iprot,ires2) = NResUreyBradleys(iprot,ires2) + 1
        if (NResUreyBradleys(iprot,ires2).gt.MAXUREYBRADLEYSPERRESIDUE) then
          write(printstring,*) "Exceeded maximum number of Urey-Bradley terms in residue ",ires2
          call printLine(WARNING,printstring)
          call printLine(WARNING,"Increase MAXUREYBRADLEYSPERRESIDUE and recompile")
          NResUreyBradleys(iprot,ires2) = NResUreyBradleys(iprot,ires2) - 1
          return
        endif
        ResidueUreyBradley(iprot,ires2,NResUreyBradleys(iprot,ires2)) = NUreyBradleys(iprot)
      endif
      if (ires1.ne.ires3 .and. ires2.ne.ires3) then
        NResUreyBradleys(iprot,ires3) = NResUreyBradleys(iprot,ires3) + 1
        if (NResUreyBradleys(iprot,ires3).gt.MAXUREYBRADLEYSPERRESIDUE) then
          write(printstring,*) "Exceeded maximum number of Urey-Bradley terms in residue ",ires3
          call printLine(WARNING,printstring)
          call printLine(WARNING,"Increase MAXUREYBRADLEYSPERRESIDUE and recompile")
          NResUreyBradleys(iprot,ires3) = NResUreyBradleys(iprot,ires3) - 1
          return
        endif
        ResidueUreyBradley(iprot,ires3,NResUreyBradleys(iprot,ires3)) = NUreyBradleys(iprot)
      endif

c     now record the residues involved in the bond
      UreyBradley(iprot,NUreyBradleys(iprot),5) = ires1
      UreyBradley(iprot,NUreyBradleys(iprot),6) = ires2
      UreyBradley(iprot,NUreyBradleys(iprot),7) = ires3
      
      return
      end
