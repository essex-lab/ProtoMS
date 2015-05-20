      subroutine uppercase(line)
      implicit none
c###################################################################
c
c     This function converts the contents of the string 'line' to
c     all uppercase
c
c     (C) Christopher Woods June 2010
c
c###################################################################

      character*(*) line
      integer lgth,i

      lgth = len(line)

      do i=1,lgth
        if (line(i:i).ge.'a' .and. line(i:i).le.'z') then
          line(i:i) = char(ichar(line(i:i))-32)
        endif
      enddo

      return
      end


      double precision function calcGeneralBond(at1, at2)
      implicit none
      double precision at1(3), at2(3)
      
      calcGeneralBond = sqrt( (at1(1)-at2(1))**2 + 
     .                        (at1(2)-at2(2))**2 +
     .                        (at1(3)-at2(3))**2 )

      return
      end

      double precision function calcGeneralAngle(at1, at2, at3)
      implicit none
      double precision at1(3), at2(3), at3(3)
      double precision vec1(3), vec2(3), id1, id2, ang
      integer i

      double precision dot, veclength2

      do i=1,3
        vec1(i) = at1(i) - at2(i)
        vec2(i) = at3(i) - at2(i)
      enddo

      id1 = vecLength2(vec1)
      id2 = vecLength2(vec2)
      
      ang = dot(vec1,vec2) / dsqrt(id1*id2)
      if (ang.gt.1.0) ang = 1.0
      if (ang.lt.-1.0) ang = -1.0
      calcGeneralAngle = dacos(ang)

      return
      end

      subroutine writeZMatLine(zid, atmnam, ipar, bndid, bnd, angid, ang,
     .                         dihid, dih, resnam, resnum)
      include 'dimensions.inc'
      include 'printer.inc'
      integer zid, ipar, bndid, angid, dihid, resnum
      double precision bnd, ang, dih, rad2Deg
      character*4 atmnam, resnam

      call uppercase(atmnam)
      call uppercase(resnam)

      write (printstring,10) zid+3, atmnam, ipar, ipar, bndid+3, bnd, angid+3,
     .                       rad2Deg(ang), dihid+3, rad2Deg(dih), resnam, resnum

10    format(i4,1x,a3,1x,i4,1x,i4,1x,i4,f12.6,i4,f12.6,i4,f12.6,1x,a3,1x,i4)

      call printLine(MCPRO,printstring)

      call lowercase(atmnam)
      call lowercase(resnam)

      return
      end


      subroutine WriteMcproBondParameters
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'printer.inc'

      integer i,j,ipar

      call prettyLine(MCPRO,"BOND PARAMETERS (COPY INTO SB FILE)")

      do i=1,NBndParam
        ipar = BndParamIndex(i)

        do j=1,2
          call uppercase(BndName(i,j))
        enddo

        write(printstring,31) BndName(i,1),BndName(i,2),BndParam(ipar,1),BndParam(ipar,2)
31      format(a2,"-",a2,1x,f5.1,4x,f8.3)
        call printLine(MCPRO,printstring)

        do j=1,2
          call uppercase(BndName(i,j))
        enddo
      enddo

      end

      subroutine WriteMcproAngleParameters
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'printer.inc'

      integer i,j,ipar
      double precision rad2Deg

      call prettyLine(MCPRO,"ANGLE PARAMETERS (COPY INTO SB FILE)")

      do i=1,NAngParam
        ipar = AngParamIndex(i)

        do j=1,3
          call uppercase(AngName(i,j))
        enddo

        write(printstring,30) AngName(i,1),AngName(i,2),AngName(i,3),
     .                        AngParam(ipar,1), rad2Deg(AngParam(ipar,2))
30      format(a2,"-",a2,"-",a2,3x,f5.1,4x,f8.2)
        call printLine(MCPRO,printstring)

        do j=1,3
          call lowercase(AngName(i,j))
        enddo
      enddo

      end

      subroutine WriteMcproDihedralParameters
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'printer.inc'

      double precision v(4)
      integer f(4)
      integer i,j,k,ipar

      call prettyLine(MCPRO,"DIHEDRAL PARAMETERS (COPY INTO PAR FILE)")

      do i=1,MAXDIHPARAM
        if (NDihTerms(i) .gt. 4) then
          write(printstring,*) "The number of dihedral terms for parameter ",i,
     .              " is greater than 4 (",NDihTerms(i),") so we can't write ",
     .              "this parameter to an mcpro file."
          call printLine(WARNING,printstring)
        else if (NDihTerms(i).gt.0) then
          do j=1,4
            f(j) = 0
            v(j) = 0
          enddo

          do j=1,NDihTerms(i)
            ipar = DihTermList(i,j)
            k = DihTerm(ipar,3)
            v(k) = 2 * DihTerm(ipar,1)
            f(k) = 0
          enddo

          write(printstring,50) i,(v(j),f(j),j=1,4)
50        format(i3,3x,f7.4,1x,i1,1x,f7.4,1x,i1,1x,f7.4,1x,i1,1x,f7.4,1x,i1)
          call printLine(MCPRO,printstring)
        endif
      enddo

      end


      subroutine WriteProteinMCPro(iprot)
      include 'dimensions.inc'
      include 'proteins.inc'  
      include 'parameters.inc'
      include 'printer.inc'    
c####################################################################
c
c     This function writes the z-matrix and parameter files
c     necessary to calculate the energy of the passed protein
c
c     (C) Christopher Woods  June 2010
c
c####################################################################

      integer iprot, i, j, k, iclj, cljID, nats
      double precision dm1(3), dm2(3), dm3(3), coords(3)
      double precision bnd, ang, dih
      double precision calcGeneralBond, calcGeneralAngle
      double precision calcGeneralDihedral, rad2Deg

      do i=1,3
        dm1(i) = 0
        dm2(i) = 0
        dm3(i) = 0
      enddo

      dm2(1) = 1
      dm3(2) = 1

      if (iprot.le.0 .or. iprot.gt.NProteins) then
        write (printstring,*) "Invalid protein number",iprot
        call printLine(WARNING,printstring)
        return
      endif

      call printLine(INFO,"Writing protein MCPRO files...")

      write (printstring,*) " ProtoMS protein z-matrix", iprot
      call printLine(MCPRO,printstring)

c     write the dummy atoms
      write (printstring,*) 
     .   "  1 DU1   -1   -1    0    0.000000   0    0.000000   0    0.000000 UNK    0"
      call printLine(MCPRO,printstring)
      write (printstring,*)
     .   "  2 DU2   -1   -1    1    1.000000   0    0.000000   0    0.000000 UNK    0"
      call printLine(MCPRO,printstring)
      write(printstring,*)
     .   "  3 DU3   -1   -1    2    1.000000   1   90.000000   0    0.000000 UNK    0"
      call printLine(MCPRO,printstring)

c     loop through all of the residues - we will build the backbone atoms
c     direct from the dummy atoms
      do i=1,NResidues(iprot)
        do j=1,3
          do k=1,3
             coords(k) = BBCoord(iprot,i,j,k)
          enddo

          bnd = calcGeneralBond(coords, dm3)
          ang = calcGeneralAngle(coords, dm3, dm2)
          dih = calcGeneralDihedral(coords, dm3, dm2, dm1)

          iclj = ZPar(iprot,BBZid(iprot,i,j))
          call writeZMatLine(BBZid(iprot,i,j), BBNam(iprot,i,j), cljID(iclj),
     .                       0, bnd, -1, ang, -2, dih, ResNam(iprot,i), ResNum(iprot,i))
        enddo

        nats = NSCAtoms(iprot,i)

        do j=1,nats
          iclj = ZPar(iprot,SCZid(iprot,i,j))
        
          call writeZMatLine(SCZid(iprot,i,j), SCNam(iprot,i,j), cljID(iclj),
     .                        ZMatSC(iprot,i,j,1), ZValSC(iprot,i,j,1),
     .                        ZMatSC(iprot,i,j,2), ZValSC(iprot,i,j,2),
     .                        ZMatSC(iprot,i,j,3), ZValSC(iprot,i,j,3),
     .                        ResNam(iprot,i), ResNum(iprot,i))
        enddo

        do j=4,4
          do k=1,3
             coords(k) = BBCoord(iprot,i,j,k)
          enddo

          bnd = calcGeneralBond(coords, dm3)
          ang = calcGeneralAngle(coords, dm3, dm2)
          dih = calcGeneralDihedral(coords, dm3, dm2, dm1)
      
          iclj = ZPar(iprot,BBZid(iprot,i,j))
          call writeZMatLine(BBZid(iprot,i,j), BBNam(iprot,i,j), cljID(iclj),
     .                       0, bnd, -1, ang, -2, dih, ResNam(iprot,i), ResNum(iprot,i))
        enddo
      enddo

      write(printstring,*) "                    Geometry Variations follow    (2I4,F12.6)"
      call printLine(MCPRO,printstring)

      write(printstring,*) "                    Variable Bonds follow         (I4)"
      call printLine(MCPRO,printstring)

      do i=1,NResidues(iprot)
        do j=1,NSCAtoms(iprot,i)
          if (ZMovSC(iprot,i,j,1)) then
            write(printstring,11) SCZid(iprot,i,j)+3
11          format(i4)
            call printLine(MCPRO,printstring)
          endif
        enddo
      enddo

      write(printstring,*) "                    Additional Bonds follow       (2I4)"
      call printLine(MCPRO,printstring)

      do i=1,NBonds(iprot)
        if (Bond(iprot,i,3).ne.0) then
          write(printstring,15) Bond(iprot,i,1)+3, Bond(iprot,i,2)+3
15        format(2i4)
          call printLine(MCPRO,printstring)
        endif
      enddo

      write(printstring,*) "                    Harmonic Constraints follow   (2I4,4F10.4)"
      call printLine(MCPRO,printstring)

      write(printstring,*) "                    Variable Bond Angles follow   (I4)"
      call printLine(MCPRO,printstring)

      do i=1,NResidues(iprot)
        do j=1,NSCAtoms(iprot,i)
          if (ZMovSC(iprot,i,j,2)) then
            write(printstring,11) SCZid(iprot,i,j)+3
            call printLine(MCPRO,printstring)
          endif
        enddo
      enddo

      write(printstring,*) "                    Additional Bond Angles follow (3I4)"
      call printLine(MCPRO,printstring)

      do i=1,NAngles(iprot)
        if (Angle(iprot,i,4).ne.0) then
          write(printstring,14) Angle(iprot,i,1)+3, Angle(iprot,i,2)+3, Angle(iprot,i,3)+3
14        format(3i4)
          call printLine(MCPRO,printstring)
        endif
      enddo

      write(printstring,*) "                    Variable Dihedrals follow     (3I4,F12.6)"
      call printLine(MCPRO,printstring)

      do i=1,NResidues(iprot)
        do j=1,NSCAtoms(iprot,i)
          if (ZMovSC(iprot,i,j,3)) then
            write(printstring,12) SCZid(iprot,i,j)+3, 0, 0, rad2Deg(ZDeltaSC(iprot,i,j,3))
12          format(3i4,f12.6)
            call printLine(MCPRO,printstring)
          endif
        enddo
      enddo

      write(printstring,*) "                    Additional Dihedrals follow   (6I4)"
      call printLine(MCPRO,printstring)

      do i=1,NDihedrals(iprot)
        if (Dihedral(iprot,i,5).ne.0) then
          write(printstring,16) (Dihedral(iprot,i,j),j=1,5),Dihedral(iprot,i,5)
16        format(6i4)
          call printLine(MCPRO,printstring)
        endif
      enddo

      write(printstring,*) "                    Domain Definitions follow     (4I4)"
      call printLine(MCPRO,printstring)

      write(printstring,*) "                    Excluded Atoms List follows   (10I4)"
      call printLine(MCPRO,printstring)

      write(printstring,*) "                    Final blank line"
      call printLine(MCPRO,printstring)

      write(printstring,*) "                                    "
      call printLine(MCPRO,printstring)

      call prettyLine(MCPRO,"CLJ PARAMETERS (COPY INTO PAR FILE)")
      call printLine(MCPRO,"TYP AN AT   CHARGE     SIGMA    EPSILON")

      do i=1,NCLJTypes
        write(printstring,20) CLJType(i), ProtonNumDB(i), AtmDB(i), ChgDB(i), SigDB(i), EpsDB(i)
20      format(i3,1x,i2,1x,a2,2x,f9.5,1x,f9.5,1x,f9.5)
        call printLine(MCPRO,printstring)
      enddo

      call writeMcproBondParameters()
      call writeMcproAngleParameters()
      call writeMcproDihedralParameters()

      end

      subroutine WriteSoluteMCPro(isol)
      include 'dimensions.inc'
      include 'printer.inc'
c####################################################################
c
c     This function writes the z-matrix and parameter files
c     for the specified solute
c
c     (C) Christopher Woods  June 2010
c
c####################################################################

      integer isol

      call printLine(WARNING, "Writing a solute MCPRO file set is not supported.")
      end

      subroutine WriteMCPro(ityp, imol)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
c####################################################################
c
c     This function writes the z-matrix and parameter files
c     to calculate the energy of the passed protein or solute
c
c     (C) Christopher Woods  June 2010
c
c####################################################################

      integer ityp, imol

      write(printstring,*) ityp, imol
      call printLine(INFO,printstring)

      if (ityp.eq.PROTEIN) then
        call writeProteinMCPro(imol)
      else
        call writeSoluteMCPro(imol)
      endif

      end
