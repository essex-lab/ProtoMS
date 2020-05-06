      double precision function calcAngle(ityp,ival,izat1,izat2,izat3)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'enums.inc'
c######################################################################
c
c     This function returns the angle between atoms izat1-izat2-izat3
c     in radians.
c
c     (C) Christopher Woods    26 January 2002
c
c######################################################################

      integer iprot,izat1,izat2,izat3,isol,ityp,ival
      double precision vec1(3),vec2(3)
      double precision ang,ia(3),ib(3),ic(3)
      double precision vecLength2,dot
      integer i
      double precision id1,id2
      
      if (ityp.eq.PROTEIN) then
        iprot = ival
        do i=1,3
          vec1(i) = ZCoords(iprot,izat1,i) - ZCoords(iprot,izat2,i)
          vec2(i) = ZCoords(iprot,izat3,i) - ZCoords(iprot,izat2,i)
        enddo
      else if (ityp.eq.SOLUTE) then
        isol = ival
        if (izat1.lt.0) then
          do i=1,3
            ia(i) = SolDummy(isol,-izat1,i)
          enddo
        else
          do i=1,3
            ia(i) = ZSolCoords(isol,izat1,i)
          enddo
        endif
        if (izat2.lt.0) then
          do i=1,3
            ib(i) = SolDummy(isol,-izat2,i)
          enddo
        else
          do i=1,3
            ib(i) = ZSolCoords(isol,izat2,i)
          enddo
        endif
        if (izat3.lt.0) then
          do i=1,3
            ic(i) = SolDummy(isol,-izat3,i)
          enddo
        else
          do i=1,3
            ic(i) = ZSolCoords(isol,izat3,i)
          enddo
        endif
        do i=1,3
          vec1(i) = ia(i) - ib(i)
          vec2(i) = ic(i) - ib(i)
        enddo
      else if (ityp.eq.ATOMS) then
        do i=1,3
          vec1(i) = TempAtoms(1,i) - TempAtoms(2,i)
          vec2(i) = TempAtoms(3,i) - TempAtoms(2,i)
        enddo
      else
        write(printstring,*) "Unknown system requested in calcAngle() - ",ityp
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is probably a programming bug")
        calcAngle = ZERO
        return
      endif          

      id1 = vecLength2(vec1)
      id2 = vecLength2(vec2)
      
      ang = dot(vec1,vec2) / dsqrt(id1*id2)
      if (ang.gt.1.0) ang = 1.0
      if (ang.lt.-1.0) ang = -1.0
      calcAngle = dacos(ang)
      
      return
      
      end
