      subroutine addMolAtom(isol,tmol,iat)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'solutes.inc'
      include 'enums.inc'
c#########################################################
c
c     This routine adds the template atom information
c     to the solute - this allows the 'zmatrix' to be
c     built, and also assigns the solute parameters
c
c     (C) Christopher Woods   19 February 2002
c
c#########################################################

      integer isol,tmol,iat,ipar1,ipar2
      character*4 iatm,ires,bndatm,bndres,angatm,angres,dihatm,dihres
      integer zatm,zbnd,zang,zdih
      integer findSoluteAtom,addCLJ,getCLJ
      double precision calcBond,calcAngle,calcDihedral

      iatm = MolTemplateAtmName(tmol,iat)
      ires = MolTemplateAtmRes(tmol,iat)
      bndatm = MolTemplateZNam(tmol,iat,1)
      bndres = MolTemplateZRes(tmol,iat,1)
      angatm = MolTemplateZNam(tmol,iat,2)
      angres = MolTemplateZRes(tmol,iat,2)
      dihatm = MolTemplateZNam(tmol,iat,3)
      dihres = MolTemplateZRes(tmol,iat,3)
      
c     first find the atoms
      zatm = findSoluteAtom(isol,iatm,ires)
      zbnd = findSoluteAtom(isol,bndatm,bndres)
      zang = findSoluteAtom(isol,angatm,angres)
      zdih = findSoluteAtom(isol,dihatm,dihres)

c     the first atom cannot be one of the special dummies
      if (zatm.lt.0) then
        write(printstring,*) "Cannot have a geometry line involving the special dummy atom ",
     .                        iatm,"-",ires
        call printLine(WARNING,printstring)
        return
      endif
      
c     all of these *must* be assigned - and decrease (as zmatricies are
c     built in sequence!)
      if (zatm.eq.0 .or. zbnd.eq.0 .or. zang.eq.0 .or. zdih.eq.0) then
        write(printstring,10) iatm,ires,zatm,bndatm,bndres,zbnd,angatm,angres,
     .                             zang,dihatm,dihres,zdih,isol,SoluteName(isol)
        call printLine(WARNING,printstring)
        call printLine(WARNING,"The missing atoms have a number 0 in brackets, e.g. (  0)")
        call printLine(WARNING,"Check to see if the missing atoms are in the coordinate file")
10      format("Problem assigning ",4(a,"-",a,"(",i4,") "),i4," ",a)
        return
      else if (zatm.le.zbnd .or. zatm.le.zang .or. zatm.le.zdih) then
        write(printstring,10) iatm,ires,zatm,bndatm,bndres,zbnd,angatm,angres,
     .                             zang,dihatm,dihres,zdih,isol,SoluteName(isol)
        call printLine(WARNING,printstring)
        call printLine(WARNING,"Check to see if all of the atoms are in the correct sequence")
        return
      endif

c     now make these the zmatrix types for this solute atom
      ZSolMat(isol,zatm,1) = zbnd
      ZSolMat(isol,zatm,2) = zang
      ZSolMat(isol,zatm,3) = zdih
      ZSolVal(isol,zatm,1) = calcBond(SOLUTE,isol,zatm,zbnd)
      ZSolVal(isol,zatm,2) = calcAngle(SOLUTE,isol,zatm,zbnd,zang)
      ZSolVal(isol,zatm,3) = calcDihedral(SOLUTE,isol,zatm,zbnd,zang,zdih)
      
c     get the clj parameter type for each atom
      ipar1 = addCLJ(MolTemplateAtmCLJ(tmol,iat,1))
      ipar2 = addCLJ(MolTemplateAtmCLJ(tmol,iat,2))
      
      ZSolPar(isol,zatm) = getCLJ(ipar1,ipar2)
      return
      end
