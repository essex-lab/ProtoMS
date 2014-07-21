      subroutine updateSoluteInternals(isol)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'solutes.inc'
c#####################################################
c
c     This routine recalculates all of the bonds/angles/
c     dihedrals for solute isol
c
c     Eventually some way must be found to cache all of this...
c
c     (C) Christopher Woods  28 February 2002
c
c#####################################################

      integer isol
      integer i,inbnds,inangs,indihs,inubys
      double precision calcBond,calcAngle,calcDihedral,calcUreyBradley
      
      if (isol.le.0 .or. isol.gt.MAXSOLUTES) return
      
c     update all of the bonds 
      inbnds = NSolBonds(isol)
      if (inbnds.gt.0) then
        do i=1,inbnds
c         do not bother if it is a zero parameter bond...
          if (SolBond(isol,i,3).ne.0) then
            SolBondLength(isol,i) = calcBond(SOLUTE,isol,SolBond(isol,i,1),SolBond(isol,i,2))
          endif
        enddo
      endif
c     now the angles
      inangs = NSolAngles(isol)
      if (inangs.gt.0) then
        do i=1,inangs
c         do not bother if it is a zero parameter angle
          if (SolAngle(isol,i,4).ne.0) then
            SolAngleSize(isol,i) = calcAngle(SOLUTE,isol,SolAngle(isol,i,1),
     .                                    SolAngle(isol,i,2),SolAngle(isol,i,3))   
          endif
        enddo
      endif
c     now the dihedrals
      indihs = NSolDihedrals(isol)
      if (indihs.gt.0) then
        do i=1,indihs
c         do not bother if it is a zero parameter dihedral
          if (SolDihedral(isol,i,5).ne.0) then
            SolDihedralSize(isol,i) = calcDihedral(SOLUTE,isol,SolDihedral(isol,i,1),
     .              SolDihedral(isol,i,2),SolDihedral(isol,i,3),SolDihedral(isol,i,4))
          endif
        enddo
      endif
c     now the Urey-Bradley terms
      inubys = NSolUreyBradleys(isol)
      if (inubys.gt.0) then
        do i=1,inubys
          if (SolUreyBradley(isol,i,3).ne.0) then
            SolUreyBradleySize(isol,i) = calcUreyBradley(SOLUTE,isol,
     .                SolUreyBradley(isol,i,1),SolUreyBradley(isol,i,2),SolUreyBradley(isol,i,3))
          endif
        enddo
      endif
      
      return
      end
