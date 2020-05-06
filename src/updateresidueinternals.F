      subroutine updateResidueInternals(iprot,ires)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'proteins.inc'
c##########################################################
c
c     This subroutine recalculates all of the bonds/angles
c     and dihedrals for residue ires in protein iprot...
c
c     Eventually some way must be found to cache this 
c     information...
c
c     (C) Christopher Woods  27 January 2002
c
c##########################################################

      integer iprot,ires
      integer i,inbnds,ibnd,inangs,iang,indihs,idih,inubys,iuby
      double precision calcBond,calcAngle,calcDihedral,calcUreyBradley

      if (iprot.le.0 .or. iprot.gt.MAXPROTEINS .or. ires.le.0) return
      
      if (ires.gt.NResidues(iprot)) return
      
c     first update all of the bonds...
      inbnds = NResBonds(iprot,ires)
      if (inbnds.gt.0) then
        do i=1,inbnds
          ibnd = ResidueBond(iprot,ires,i)
          BondLength(iprot,ibnd) = calcBond(PROTEIN,iprot,Bond(iprot,ibnd,1),Bond(iprot,ibnd,2))
        enddo
      endif
      
c     now update all of the angles
      inangs = NResAngles(iprot,ires)
      if (inangs.gt.0) then
        do i=1,inangs
          iang = ResidueAngle(iprot,ires,i)
          AngleSize(iprot,iang) = calcAngle(PROTEIN,iprot,Angle(iprot,iang,1),Angle(iprot,iang,2),
     .                                                  Angle(iprot,iang,3))
        enddo
      endif
      
c     now update all of the dihedrals
      indihs = NResDihedrals(iprot,ires)
      if (indihs.gt.0) then
        do i=1,indihs
          idih = ResidueDihedral(iprot,ires,i)
          DihedralSize(iprot,idih) = calcDihedral(PROTEIN,iprot,Dihedral(iprot,idih,1),Dihedral(iprot,idih,2),
     .                                                     Dihedral(iprot,idih,3),Dihedral(iprot,idih,4))
        enddo
      endif
      
c     now update all of the Urey-Bradley terms
      inubys = NResUreyBradleys(iprot,ires)
      if (inubys.gt.0) then
        do i=1,inubys
          iuby = ResidueUreyBradley(iprot,ires,i)
          UreyBradleySize(iprot,iuby) = calcUreyBradley(PROTEIN,iprot,
     .                                UreyBradley(iprot,iuby,1),UreyBradley(iprot,iuby,2))
        enddo
      endif
      
      return 
      
      end
