      subroutine addProteinBond(iprot,iat1,iat2,dummy,flex,delta,ipar0,ipar1)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'proteins.inc'
c######################################################
c
c     This function adds a bond between atoms with zmatrix
c     id's of zat1 and zat2...
c
c     if typ = 0 then it is a dummy bond...
c     if flx = 1 then it is a flexible bond...
c           with a range of motion of +/- delta
c
c     (C) Christopher Woods   30 January 2002
c
c######################################################

      integer iprot,zat1,zat2,ires1,ires2,nbnds,iparam
      integer iat1,iat2,ipar0,ipar1
      logical dummy,flex
      integer getBond,setDirectBond
      double precision delta
      logical hasBond
      
      zat1 = iat1
      zat2 = iat2
c     sort the atoms
      call sortZBond(zat1,zat2)
      
      ires1 = Zid2ResID(iprot,zat1,1)
      ires2 = Zid2ResID(iprot,zat2,1)

c     make sure that the bond hasn't already been added...
      if (hasBond(iprot,zat1,zat2)) then
        return
      endif
      
c     now add the bond...
      NBonds(iprot) = NBonds(iprot) + 1
      nbnds = NBonds(iprot)
      if (nbnds.gt.MAXBONDS) then
        write(printstring,*) "Exceeded maximum number of protein bonds! (",MAXBONDS,")"
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Increase MAXBONDS and recompile")
        NBonds(iprot) = NBonds(iprot) - 1
        return
      endif
      Bond(iprot,nbnds,1) = zat1
      Bond(iprot,nbnds,2) = zat2

      BondDummy(iprot,NBonds(iprot)) = .false.
      if (flex) then
        FlexBond(iprot,nbnds) = .true.
        call makeVariableBond(iprot,zat1,zat2,delta)
      endif
      
      if (dummy) then
c       this bond is a dummy and has no energy parameters
        BondDummy(iprot,NBonds(iprot)) = .true.
        iparam = 0
      else if (ipar0.ge.0 .and. ipar1.ge.0) then
        iparam = setDirectBond(ipar0,ipar1)
      else
        iparam = getBond(ZPar(iprot,zat1),ZPar(iprot,zat2))
      endif
      Bond(iprot,NBonds(iprot),3) = iparam

c     now add the bond to the residue
      NResBonds(iprot,ires1) = NResBonds(iprot,ires1) + 1
      if (NResBonds(iprot,ires1).gt.MAXBONDSPERRESIDUE) then
        write(printstring,*) "Exceeded maximum number of bonds in residue ",ires1
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Increase MAXBONDSPERRESIDUE and recompile")
        NResBonds(iprot,ires1) = NResBonds(iprot,ires1) - 1
        return
      endif
      ResidueBond(iprot,ires1,NResBonds(iprot,ires1)) = NBonds(iprot)
      if (ires1.ne.ires2) then
        NResBonds(iprot,ires2) = NResBonds(iprot,ires2) + 1
        ResidueBond(iprot,ires2,NResBonds(iprot,ires2)) = NBonds(iprot)
        if (NResBonds(iprot,ires2).gt.MAXBONDSPERRESIDUE) then
          write(printstring,*) "Exceeded maximum number of bonds in residue ",ires2
          call printLine(WARNING,printstring)
          call printLine(WARNING,"Increase MAXBONDSPERRESIDUE and recompile")
          NResBonds(iprot,ires1) = NResBonds(iprot,ires1) - 1
          return
        endif
      endif

c     now record the residues involved in the bond
      Bond(iprot,NBonds(iprot),4) = ires1
      Bond(iprot,NBonds(iprot),5) = ires2
      
      return
      end
