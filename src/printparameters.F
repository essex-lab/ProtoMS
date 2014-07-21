
      character*4 function SolAtomName(isol, iatm)
      include 'dimensions.inc'
      include 'solutes.inc'

      integer isol, iatm

      if (iatm.gt.0) then
          SolAtomName = SolAtomNam(isol,iatm)
      else if (iatm.eq.-1) then
          SolAtomName = "DM1"
      else if (iatm.eq.-2) then
          SolAtomName = "DM2"
      else if (iatm.eq.-3) then
          SolAtomName = "DM3"
      else
          SolAtomName = "DUM"
      endif

      return
      end

      character*4 function SolResidueName(isol, iatm)
      include 'dimensions.inc'
      include 'solutes.inc'

      integer isol, iatm

      if (iatm.gt.0) then
        SolResidueName = SolResNam( isol, SolResidue(isol,iatm) )
      else
        SolResidueName = "DUM"
      endif
 
      return
      end

      character*8 function SolPertType(typ)
      implicit none

      integer typ

      if (typ .eq. 1) then
        SolPertType = "BOND    "
      else if (typ .eq. 2) then
        SolPertType = "ANGLE   "
      else
        SolPertType = "DIHEDRAL"
      endif

      return
      end

      subroutine PrintParameters()
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'templates.inc'
      include 'simparams.inc'
      include 'parameters.inc'

c     #######################################
c     #
c     # This subroutine prints out all of the
c     # parameters used for all of the loaded
c     # molecules to a Parameters stream. You
c     # can use this to check to see which parameters
c     # are being used, or as a means by which you
c     # can use ProtoMS to parameterise molecules
c     # (e.g. this is used by the ProtoMS module
c     # in Sire)
c     #
c     # (C) Christopher Woods, November 2008
c     #
c     #######################################

      double precision Charge, Sigma, Epslon
      double precision ir,ik
      integer zp

      character*4 SolAtomName, SolResidueName
      character*4 AtomName
      character*8 SolPertType

      double precision rad2Deg

      integer getResResBondIndex

      integer iat1,iat2,iat3,iat4,ipar,idx
      integer i,j,k,l

      call prettyLine(PPARAMS, "Writing out molecular parameters")

      do i=1,NProteins
        write(printstring,*) "Protein ",i," ",ProteinName(i)
        call printLine(PPARAMS, printstring)

c       print out the zmatrix
        do j=1,NResidues(i)
          do k=1,NSCAtoms(i,j)        
            iat1 = SCZid(i,j,k)
            iat2 = ZMatSC(i,j,k,1)
            iat3 = ZMatSC(i,j,k,2)
            iat4 = ZMatSC(i,j,k,3)

            write(printstring,100) AtomName(i,iat1), ResNum(i,Zid2ResID(i,iat1,1)), 
     .                             AtomName(i,iat2), ResNum(i,Zid2ResID(i,iat2,1)),
     .                             AtomName(i,iat3), ResNum(i,Zid2ResID(i,iat3,1)), 
     .                             AtomName(i,iat4), ResNum(i,Zid2ResID(i,iat4,1))

            call printLine(PPARAMS, printstring)
 
100         format("ZMATRIX ",a4," ( ",i5," ) - ",
     .                        a4," ( ",i5," ) - ",
     .                        a4," ( ",i5," ) - ",
     .                        a4," ( ",i5," )")

          enddo
        enddo

c       print out the CLJ parameters
        do j=1,NZAtoms(i)
          write(printstring,101), AtomName(i,j), ResNum(i,Zid2ResID(i,j,1)),
     .                            Charge(ZPar(i,j)), Sigma(ZPar(i,j)),
     .                            Epslon(ZPar(i,j))
          call printLine(PPARAMS, printstring)
          
101       format("Atom ",a4," ( ",i5," ) Charge= ",f12.7," Sigma= ",f12.7," Epsilon= ",f12.7)

          if (ZPar(i,j) .lt. 0) then
            zp = -ZPar(i,j)
            write(printstring,1001) AtomName(i,j), ResNum(i,Zid2ResID(i,j,1)),
     .                              VariableChg(zp,1), VariableChg(zp,3),
     .                              VariableSig(zp,1), VariableSig(zp,3),
     .                              VariableEps(zp,1), VariableEps(zp,3)
            call printLine(PPARAMS, printstring)

1001        format("AtomPert ",a4," ( ",i5," ) Charge= ",f12.7," ",f12.7,
     .             " Sigma= ",f12.7," ",f12.7," Epsilon= ",f12.7," ",f12.7)

          endif
        enddo

c       print out the bonded pairs and scale factors (non-bonded pairs would be too much!)
        do j=1, NResidues(i)
          do k=1,NResNB(i,j)
            iat1 = ResNB(i,j,k,1)
            iat2 = ResNB(i,j,k,2)
            write(printstring,125) AtomName(i,iat1), ResNum(i,Zid2ResID(i,iat1,1)),
     .                             AtomName(i,iat2), ResNum(i,Zid2ResID(i,iat2,1)),
     .                             1.0, 1.0
            call printLine(PPARAMS,printstring)
          enddo 

          do k=1,NResNB14(i,j)
            iat1 = ResNB14(i,j,k,1)
            iat2 = ResNB14(i,j,k,2)
            write(printstring,125) AtomName(i,iat1), ResNum(i,Zid2ResID(i,iat1,1)),
     .                             AtomName(i,iat2), ResNum(i,Zid2ResID(i,iat2,1)),
     .                             OneOverScl14C, OneOverScl14LJ
            call printLine(PPARAMS,printstring)
          enddo
        enddo

c       print out the bonded pairs and scale factors for the residue-residue interactions
        do j=1, NResidues(i) - 1
          do k=j+1, NResidues(i)

            idx = getResResBondIndex(i,j,k)

            if (idx .gt. 0) then
              do l=1,NResResNB(i,idx)
                iat1 = ResResNB(i,idx,l,1)
                iat2 = ResResNB(i,idx,l,2)
                write(printstring,125) AtomName(i,iat1), ResNum(i,Zid2ResID(i,iat1,1)),
     .                                 AtomName(i,iat2), ResNum(i,Zid2ResID(i,iat2,1)),
     .                                 0.0, 0.0
                call printLine(PPARAMS,printstring)
              enddo

              do l=1,NResResNB14(i,idx)
                iat1 = ResResNB14(i,idx,l,1)
                iat2 = ResResNB14(i,idx,l,2)
                 write(printstring,125) AtomName(i,iat1), ResNum(i,Zid2ResID(i,iat1,1)),
     .                                  AtomName(i,iat2), ResNum(i,Zid2ResID(i,iat2,1)),
     .                                  OneOverScl14C, OneOverScl14LJ
                 call printLine(PPARAMS,printstring)
              enddo
            endif
          enddo
        enddo

125     format("NB ", a4," ( ",i4," ) - ",a4," ( ",i4," ) CSCL= ",f12.7," LJSCL= ",f12.7)

        write(printstring,*) "nBond, angle, dihedral ",NBonds(i)," ",NAngles(i)," ",NDihedrals(i)
        call printLine(INFO,printstring)

c       print out the bond parameters
        do j=1, NBonds(i)
          iat1 = Bond(i,j,1)
          iat2 = Bond(i,j,2)
          ipar = Bond(i,j,3)

          if (ipar.lt.0) then
            ik = VariableBndParam(-ipar,2,1)
            ir = VariableBndParam(-ipar,2,2)
          else if (ipar.eq.0) then
            ik = 0
            ir = 0
          else
            ik = BndParam(ipar,1)
            ir = BndParam(ipar,2)
          endif

          write(printstring,103), AtomName(i,iat1), ResNum(i,Zid2ResID(i,iat1,1)),
     .                            AtomName(i,iat2), ResNum(i,Zid2ResID(i,iat2,1)),
     .                            ik, ir
          call printLine(PPARAMS, printstring)

          if (.not. BondDummy(i,j)) then
            write(printstring,1103), AtomName(i,iat1), ResNum(i,Zid2ResID(i,iat1,1)),
     .                               AtomName(i,iat2), ResNum(i,Zid2ResID(i,iat2,1))
            call printLine(PPARAMS, printstring)
          endif

103       format("Bond ",a4," ( ",i5," ) - ",a4," ( ",i5," ) K= ",f12.7," R= ",f12.7)
1103      format("Connect ",a4," ( ",i5," ) - ",a4," ( ",i5," )")
        enddo

c       print out the bond deltas
        do j=1,NResidues(i)
          do k=1,NSCAtoms(i,j)        
            if (ZMovSC(i,j,k,1)) then
              iat1 = SCZid(i,j,k)
              iat2 = ZMatSC(i,j,k,1)
            
              write(printstring,110) AtomName(i,iat1), ResNum(i,Zid2ResID(i,iat1,1)),
     .                               AtomName(i,iat2), ResNum(i,Zid2ResID(i,iat2,1)),
     .                               ZDeltaSC(i,j,k,1)
              call printLine(PPARAMS,printstring)
            endif

110         format("BondDelta ",a4," ( ",i5," ) - ",a4," ( ",i5," ) Delta= ",f12.7)
          enddo
        enddo

c       print out the angle parameters
        do j=1,NAngles(i)
          iat1 = Angle(i,j,1)
          iat2 = Angle(i,j,2)
          iat3 = Angle(i,j,3)
          ipar = Angle(i,j,4)

          if (ipar.lt.0) then
            ik = VariableAngParam(-ipar,2,1)
            ir = VariableAngParam(-ipar,2,2)
          else if (ipar.eq.0) then
            ik = 0
            ir = 0
          else
            ik = AngParam(ipar,1)
            ir = AngParam(ipar,2)
          endif

          write(printstring,104), AtomName(i,iat1), ResNum(i,Zid2ResID(i,iat1,1)),
     .                            AtomName(i,iat2), ResNum(i,Zid2ResID(i,iat2,1)),
     .                            AtomName(i,iat3), ResNum(i,Zid2ResID(i,iat3,1)),
     .                            ik, rad2Deg(ir)
 
          call printLine(PPARAMS, printstring)
  
104       format("Angle ",a4," ( ",i5," ) - ",a4," ( ",i5," ) - ",a4," ( ",i5," ) K= ",
     .           f12.7," R= ",f12.7)
        enddo

c       print out the angle deltas
        do j=1,NResidues(i)
          do k=1,NSCAtoms(i,j)        
            if (ZMovSC(i,j,k,2)) then
              iat1 = SCZid(i,j,k)
              iat2 = ZMatSC(i,j,k,1)
              iat3 = ZMatSC(i,j,k,2)            

              write(printstring,112) AtomName(i,iat1), ResNum(i,Zid2ResID(i,iat1,1)),
     .                               AtomName(i,iat2), ResNum(i,Zid2ResID(i,iat2,1)),
     .                               AtomName(i,iat3), ResNum(i,Zid2ResID(i,iat3,1)),
     .                               rad2Deg(ZDeltaSC(i,j,k,2))
              call printLine(PPARAMS,printstring)
            endif

112         format("AngleDelta ",a4," ( ",i5," ) - ",a4," ( ",i5," ) - ",
     .                           a4," ( ",i5," ) Delta= ",f12.7)
          enddo
        enddo

c       print out the Urey-Bradley parameters
        do j=1,NUreyBradleys(i)
          iat1 = UreyBradley(i,j,1)
          iat2 = UreyBradley(i,j,2)
          ipar = UreyBradley(i,j,3)

          if (ipar.lt.0) then
            ik = VariableUByParam(-ipar,2,1)
            ir = VariableUByParam(-ipar,2,2)
          else if (ipar.eq.0) then
            ik = 0
            ir = 0
          else
            ik = UByParam(ipar,1)
            ir = UByParam(ipar,2)
          endif

          write(printstring,106), AtomName(i,iat1), ResNum(i,Zid2ResID(i,iat1,1)),
     .                            AtomName(i,iat2), ResNum(i,Zid2ResID(i,iat2,1)),
     .                            ik, ir
          call printLine(PPARAMS, printstring)
106       format("UreyBradley ",a4," ( ",i5," ) - ",a4," ( ",i5," ) R=",
     .                            f12.7," K= ",f12.7)
        enddo

c       print out the dihedral parameters
        do j=1,NDihedrals(i)
          iat1 = Dihedral(i,j,1)
          iat2 = Dihedral(i,j,2)
          iat3 = Dihedral(i,j,3)
          iat4 = Dihedral(i,j,4)
          ipar = Dihedral(i,j,5)

          write(printstring,107), AtomName(i,iat1), ResNum(i,Zid2ResID(i,iat1,1)),
     .                            AtomName(i,iat2), ResNum(i,Zid2ResID(i,iat2,1)),
     .                            AtomName(i,iat3), ResNum(i,Zid2ResID(i,iat3,1)),
     .                            AtomName(i,iat4), ResNum(i,Zid2ResID(i,iat4,1))
          call printLine(PPARAMS, printstring)

          call printDihedralParameter(ipar)
        
107       format("Dihedral ",a4," ( ",i5," ) - ",a4," ( ",i5," ) - ",a4," ( ",i5," ) - ",
     .                       a4," ( ",i5," )")
        enddo

c       print out the dihedral deltas
        do j=1,NResidues(i)
          do k=1,NSCAtoms(i,j)        
            if (ZMovSC(i,j,k,3)) then
              iat1 = SCZid(i,j,k)
              iat2 = ZMatSC(i,j,k,1)
              iat3 = ZMatSC(i,j,k,2)            
              iat4 = ZMatSC(i,j,k,3)

              write(printstring,114) AtomName(i,iat1), ResNum(i,Zid2ResID(i,iat1,1)),
     .                               AtomName(i,iat2), ResNum(i,Zid2ResID(i,iat2,1)),
     .                               AtomName(i,iat3), ResNum(i,Zid2ResID(i,iat3,1)),
     .                               AtomName(i,iat4), ResNum(i,Zid2ResID(i,iat4,1)),
     .                               rad2Deg(ZDeltaSC(i,j,k,3))
              call printLine(PPARAMS,printstring)
            endif

114         format("DihedralDelta ",a4," ( ",i5," ) - ",a4," ( ",i5," ) - ",
     .                              a4," ( ",i5," ) - ",a4," ( ",i5," ) Delta= ",f12.7)
          enddo
        enddo

      enddo

      do i=1,NSolutes
        if (IgnoreSolute(i)) goto 1044

        write(printstring,*) "Solute ",i," ",SoluteName(i)
        call printLine(PPARAMS, printstring)

c       print out the zmatrix
        do j=1,NSolAtoms(i)
          iat1 = j
          iat2 = ZSolMat(i,j,1)
          iat3 = ZSolMat(i,j,2)
          iat4 = ZSolMat(i,j,3)
         
          write(printstring,19) SolAtomName(i,iat1), SolResidueName(i,iat1),
     .                          SolAtomName(i,iat2), SolResidueName(i,iat2),
     .                          SolAtomName(i,iat3), SolResidueName(i,iat3),
     .                          SolAtomName(i,iat4), SolResidueName(i,iat4)

          call printLine(PPARAMS,printstring)

19        format("ZMATRIX ",a4," ( ",a4," ) - ",
     .                     a4," ( ",a4," ) - ",
     .                     a4," ( ",a4," ) - ",
     .                     a4," ( ",a4," )")
        enddo

c       print out any variable geometry
        do j=1,NMolTemplVar(i)

          if (MolTemplVarTyp(i,j) .eq. 1) then
            write(printstring,1019) MolTemplVarNam(i,j), MolTemplVarRes(i,j),
     .                              SolPertType(MolTemplVarTyp(i,j)), 
     .                              MolTemplVarVal(i,j,1), MolTemplVarVal(i,j,2)
          else
            write(printstring,1019) MolTemplVarNam(i,j), MolTemplVarRes(i,j),
     .                              SolPertType(MolTemplVarTyp(i,j)),
     .                              rad2Deg(MolTemplVarVal(i,j,1)), rad2Deg(MolTemplVarVal(i,j,2))
          endif

          call printLine(PPARAMS,printstring)

1019      format("ZMATRIXPERT ",a4," ( ",a4," ) ",a8," ",f12.5," ",f12.5)

        enddo

c       print out the CLJ parameters
        do j=1,NSolAtoms(i)
          write(printstring,20), SolAtomName(i,j), SolResidueName(i,j),
     .                            Charge(ZSolPar(i,j)), Sigma(ZSolPar(i,j)),
     .                                                  Epslon(ZSolPar(i,j))
          call printLine(PPARAMS, printstring)

20        format("Atom ",a4," ( ",a4," ) Charge= ",f12.7," Sigma= ",f12.7," Epsilon= ",f12.7)

          if (ZSolPar(i,j) .lt. 0) then
            zp = -ZSolPar(i,j)

            write(printstring,1020) SolAtomName(i,j), SolResidueName(i,j),
     .                              VariableChg(zp,1), VariableChg(zp,3),
     .                              VariableSig(zp,1), VariableSig(zp,3),
     .                              VariableEps(zp,1), VariableEps(zp,3)
            call printLine(PPARAMS, printstring)

1020        format("AtomPert ",a4," ( ",a4," ) Charge= ",f12.7," ",f12.7,
     .             " Sigma=",f12.7," ",f12.7," Epsilon=",f12.7," ",f12.7)
         endif

        enddo

c       print out the bond parameters
        do j=1,NSolBonds(i)
          iat1 = SolBond(i,j,1)
          iat2 = SolBond(i,j,2)
          ipar = SolBond(i,j,3)

          if (ipar.lt.0) then
            ik = VariableBndParam(-ipar,2,1)
            ir = VariableBndParam(-ipar,2,2)
          else if (ipar.eq.0) then
            ik = 0
            ir = 0
          else
            ik = BndParam(ipar,1)
            ir = BndParam(ipar,2)
          endif

          write(printstring,21), SolAtomName(i,iat1), SolResidueName(i,iat1),
     .                           SolAtomName(i,iat2), SolResidueName(i,iat2),
     .                           ik, ir                         
          call printLine(PPARAMS, printstring)

          if (ipar.lt.0) then
            write(printstring,2021), SolAtomName(i,iat1), SolResidueName(i,iat1),
     .                               SolAtomName(i,iat2), SolResidueName(i,iat2),
     .                               VariableBndParam(-ipar,1,1), VariableBndParam(-ipar,1,2),
     .                               VariableBndParam(-ipar,3,1), VariableBndParam(-ipar,3,2)
            call printLine(PPARAMS,printstring)
2021        format("BondPert ",a4," ( ",a4," ) - ",a4," ( ",a4," ) KB= ",f12.7," RB= ",f12.8,
     .              " KF= ",f12.7," RF= ",f12.7)
         endif

          if (.not. SolBondDummy(i,j)) then
            write(printstring,1021) SolAtomName(i,iat1), SolResidueName(i,iat1),
     .                              SolAtomName(i,iat2), SolResidueName(i,iat2)
            call printLine(PPARAMS,printstring)
          endif

21        format("Bond ",a4," ( ",a4," ) - ",a4," ( ",a4," ) K= ",f12.7," R= ",f12.7)
1021      format("Connect ",a4," ( ",a4,") - ",a4," ( ",a4," )")
        enddo

c       print out the non-bonded pair scale factors
        do j=1,NSolNB(i)
          iat1 = SolNB(i,j,1)
          iat2 = SolNB(i,j,2)

          if (SolNB14(i,j)) then
            write(printstring,25) SolAtomName(i,iat1), SolResidueName(i,iat1), 
     .                            SolAtomName(i,iat2), SolResidueName(i,iat2), 
     .                            OneOverScl14C, OneOverScl14LJ

          else
            write(printstring,25) SolAtomName(i,iat1), SolResidueName(i,iat1),
     .                            SolAtomName(i,iat2), SolResidueName(i,iat2),
     .                            1.0, 1.0

          endif

          call printLine(PPARAMS,printstring)

25        format("NB ", a4," ( ",a4," ) - ",a4," ( ",a4," ) CSCL= ",f12.7," LJSCL= ",f12.7)

        enddo

c       print out the bond deltas
        do j=1,NSolAtoms(i)
          if (ZSolMov(i,j,1)) then
            iat1 = j
            iat2 = ZSolMat(i,j,1)
            
            write(printstring,41) SolAtomName(i,iat1), SolResidueName(i,iat1),
     .                            SolAtomName(i,iat2), SolResidueName(i,iat2),
     .                            ZSolDelta(i,j,1)
            call printLine(PPARAMS,printstring)
          endif

41        format("BondDelta ",a4," ( ",a4," ) - ",a4," ( ",a4," ) Delta= ",f12.7)
        enddo

c       print out the angle parameters
        do j=1,NSolAngles(i)
          iat1 = SolAngle(i,j,1)
          iat2 = SolAngle(i,j,2)
          iat3 = SolAngle(i,j,3)
          ipar = SolAngle(i,j,4)

          if (ipar.lt.0) then
            ik = VariableAngParam(-ipar,2,1)
            ir = VariableAngParam(-ipar,2,2)
          else if (ipar.eq.0) then
            ik = 0
            ir = 0
          else
            ik = AngParam(ipar,1)
            ir = AngParam(ipar,2)
          endif

          write(printstring,22)  SolAtomName(i,iat1), SolResidueName(i,iat1),
     .                           SolAtomName(i,iat2), SolResidueName(i,iat2),
     .                           SolAtomName(i,iat3), SolResidueName(i,iat3),
     .                           ik, rad2Deg(ir)
          call printLine(PPARAMS, printstring)

          if (ipar.lt.0) then
            write(printstring,2022)  SolAtomName(i,iat1), SolResidueName(i,iat1),
     .                               SolAtomName(i,iat2), SolResidueName(i,iat2),
     .                               SolAtomName(i,iat3), SolResidueName(i,iat3),
     .                               VariableAngParam(-ipar,1,1), rad2Deg(VariableAngParam(-ipar,1,2)),
     .                               VariableAngParam(-ipar,3,1), rad2Deg(VariableAngParam(-ipar,3,2))
            call printLine(PPARAMS,printstring)
2022        format("AnglePert ",a4," ( ",a4," ) - ",a4," ( ",a4," ) - ",a4," ( ",a4," ) KB= ",
     .              f12.7," RB= ",f12.8," KF= ",f12.7," RF= ",f12.7)
         endif

22        format("Angle ",a4," ( ",a4," ) - ",a4," ( ",a4," ) - ",a4," ( ",a4," ) K= ",
     .           f12.7," R= ",f12.7)
        enddo

c       print out the angle deltas
        do j=1,NSolAtoms(i)
          if (ZSolMov(i,j,2)) then
            iat1 = j
            iat2 = ZSolMat(i,j,1)
            iat3 = ZSolMat(i,j,2)
            
            write(printstring,42) SolAtomName(i,iat1), SolResidueName(i,iat1),
     .                            SolAtomName(i,iat2), SolResidueName(i,iat2),
     .                            SolAtomName(i,iat3), SolResidueName(i,iat3),
     .                            rad2Deg(ZSolDelta(i,j,2))
            call printLine(PPARAMS,printstring)
          endif
     
42        format("AngleDelta ",a4," ( ",a4," ) - ",a4," ( ",a4," ) - ",a4," ( ",a4," ) Delta= ",f12.7)
        enddo

c       print out the Urey-Bradley parameters
        do j=1,NSolUreyBradleys(i)
          iat1 = SolUreyBradley(i,j,1)
          iat2 = SolUreyBradley(i,j,2)
          ipar = SolUreyBradley(i,j,3)

          if (ipar.lt.0) then
            ik = VariableUByParam(-ipar,2,1)
            ir = VariableUByParam(-ipar,2,2)
          else if (ipar.eq.0) then
            ik = 0
            ir = 0
          else
            ik = UByParam(ipar,1)
            ir = UByParam(ipar,2)
          endif

          write(printstring,23), SolAtomName(i,iat1), SolResidueName(i,iat1),
     .                           SolAtomName(i,iat2), SolResidueName(i,iat2),
     .                           ik, ir
          call printLine(PPARAMS, printstring)

          if (ipar.lt.0) then
            write(printstring,2023), SolAtomName(i,iat1), SolResidueName(i,iat1),
     .                               SolAtomName(i,iat2), SolResidueName(i,iat2),
     .                               VariableUByParam(-ipar,1,1), VariableUByParam(-ipar,1,2),
     .                               VariableUByParam(-ipar,3,1), VariableUByParam(-ipar,3,2)
            call printLine(PPARAMS,printstring)
2023        format("UreyBradleyPert ",a4," ( ",a4," ) - ",a4," ( ",a4," ) KB= ",f12.7," RB= ",f12.8,
     .              " KF= ",f12.7," RF= ",f12.7)
         endif

23        format("UreyBradley ",a4," ( ",a4," ) - ",a4," ( ",a4," ) R=",
     .                              f12.7," K= ",f12.7)
        enddo

c       print out the dihedral parameters
        do j=1,NSolDihedrals(i)
          iat1 = SolDihedral(i,j,1)
          iat2 = SolDihedral(i,j,2)
          iat3 = SolDihedral(i,j,3)
          iat4 = SolDihedral(i,j,4)
          ipar = SolDihedral(i,j,5)

          write(printstring,24), SolAtomName(i,iat1), SolResidueName(i,iat1),
     .                           SolAtomName(i,iat2), SolResidueName(i,iat2),
     .                           SolAtomName(i,iat3), SolResidueName(i,iat3),
     .                           SolAtomName(i,iat4), SolResidueName(i,iat4)
          call printLine(PPARAMS, printstring)

          call printDihedralParameter(ipar)
        
24        format("Dihedral ",a4," ( ",a4," ) - ",a4," ( ",a4," ) - ",a4," ( ",a4," ) - ",
     .                       a4," ( ",a4," )")
        enddo

c       print out the dihedral deltas
        do j=1,NSolAtoms(i)
          if (ZSolMov(i,j,3)) then
            iat1 = j
            iat2 = ZSolMat(i,j,1)
            iat3 = ZSolMat(i,j,2)
            iat4 = ZSolMat(i,j,3)

            write(printstring,44) SolAtomName(i,iat1), SolResidueName(i,iat1),
     .                            SolAtomName(i,iat2), SolResidueName(i,iat2),
     .                            SolAtomName(i,iat3), SolResidueName(i,iat3),
     .                            SolAtomName(i,iat4), SolResidueName(i,iat4),
     .                            rad2Deg(ZSolDelta(i,j,3))
            call printLine(PPARAMS,printstring)
          endif

44        format("DihedralDelta ",a4," ( ",a4," ) - ",a4," ( ",a4," ) - ",a4," ( ",a4," ) - ",
     .                            a4," ( ",a4," ) Delta= ",f12.7)
        enddo

1044    continue
      enddo

      do i=1,NSolvents
        write(printstring,*) "Solvent ",i," ",SvnName(i)
        call printLine(PPARAMS, printstring)

c       print out the CLJ parameters
        do j=1,NSvnAtoms(i)
          write(printstring,20) SvnAtomName(i,j), SvnName(i), Charge(ZSvnPar(i,j)), 
     .                              Sigma(ZSvnPar(i,j)), Epslon(ZSvnPar(i,j))        
          call printLine(PPARAMS, printstring)

          if (ZSvnPar(i,j) .lt. 0) then
            zp = -ZSvnPar(i,j)
            write(printstring,1020), SvnAtomName(i,j), SvnName(i),
     .                               VariableChg(zp,1), VariableChg(zp,3),
     .                               VariableSig(zp,1), VariableSig(zp,3),
     .                               VariableEps(zp,1), VariableEps(zp,3)

            call printLine(PPARAMS, printstring)    

          endif
        enddo

      enddo

      return

      end

      subroutine printDihedralParameter(ipar)
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'printer.inc'      

      integer ipar,iparf,iparb,i,j

      write(printstring,*) "HELP ",ipar
      call printLine(INFO,printstring)

      if (ipar.eq.0) then
        write(printstring,10) 0.0,0.0,0.0,0.0
        call printLine(PPARAMS,printstring)
        return
      else
        if (ipar.lt.0) then
          iparb = VariableDihParam(-ipar,1)
          iparf = VariableDihParam(-ipar,2)
          ipar = iparb
        else 
          iparb = -1
          iparf = -1
        endif

        write (printstring,*) "HELP2 ",iparb," ",iparf
        call printLine(INFO,printstring)

        do i=1,NDihTerms(ipar)
          write(printstring,10) (DihTerm(DihTermList(ipar,i),j), j=1,4)  
          call printLine(PPARAMS,printstring)

10        format("DihedralParameter ik1= ",f12.7," ik2= ",f12.7," ik3= ",f12.7," ik4= ",f12.7)
        enddo

        if (iparf.eq.0) then
          write(printstring,20) 0.0,0.0,0.0,0.0
          call printLine(PPARAMS,printstring)
        else if (iparf.gt.0) then
          do i=1,NDihTerms(iparf)
            write(printstring,20) (DihTerm(DihTermList(iparf,i),j), j=1,4)
            call printLine(PPARAMS,printstring)
          
20          format("DihedralParameterPertF ik1= ",f12.7," ik2= ",f12.7," ik3= ",f12.7," ik4= ",f12.7)
          enddo
        endif

        if (iparb.eq.0) then
          write(printstring,30) 0.0,0.0,0.0,0.0
          call printLine(PPARAMS,printstring)
        else if (iparb.gt.0) then
          do i=1,NDihTerms(iparb)
            write(printstring,30) (DihTerm(DihTermList(iparb,i),j), j=1,4)
            call printLine(PPARAMS,printstring)
        
30          format("DihedralParameterPertB ik1= ",f12.7," ik2= ",f12.7," ik3= ",f12.7," ik4= ",f12.7)
          enddo
        endif

      endif

      return
      end

