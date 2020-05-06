      double precision function calcBond(ityp,ival,izat1,izat2)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'enums.inc'
c#############################################################
c
c     Ths subroutine calculates the length of the bond between
c     atoms with zmatrix id izat1 and izat2, and returns it...
c
c     (C) Christopher Woods   26 January 2002
c
c#############################################################

      integer iprot,izat1,izat2,ityp,isol,ival,i
      double precision lngth,ia(3),ib(3)

      if (ityp.eq.PROTEIN) then
        iprot = ival
        lngth = (ZCoords(iprot,izat1,1)-ZCoords(iprot,izat2,1))**2
     .          + (ZCoords(iprot,izat1,2)-ZCoords(iprot,izat2,2))**2
     .            + (ZCoords(iprot,izat1,3)-ZCoords(iprot,izat2,3))**2
      else if (ityp.eq.SOLUTE) then
        isol = ival
        if (izat1.lt.0) then
          do i=1,3
            ia(i) = SolDummy(isol,-izat1,i)
          enddo
        else if (izat1.eq.0) then
          do i=1,3
            ia(i) = ZERO
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
        else if (izat2.eq.0) then
          do i=1,3
            ib(i) = ZERO
          enddo
        else
          do i=1,3
            ib(i) = ZSolCoords(isol,izat2,i)
          enddo
        endif
        lngth = (ia(1)-ib(1))**2
     .          + (ia(2)-ib(2))**2
     .            + (ia(3)-ib(3))**2
      else if (ityp.eq.ATOMS) then
        lngth = (TempAtoms(1,1)-TempAtoms(2,1))**2
     .         +(TempAtoms(1,2)-TempAtoms(2,2))**2
     .         +(TempAtoms(1,3)-TempAtoms(2,3))**2
      else
        write(printstring,*) "Unrecognised system requested in calcBond() - ",ityp
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is probably a programming bug!")
        calcBond = ZERO
        return
      endif
     
      calcBond = sqrt(lngth)

      return
      
      end
