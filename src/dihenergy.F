      double precision function dihEnergy(ipar,ang)
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'constants.inc'
c##############################################################
c
c     This function returns the dihedral energy for parameter
c     ipar at angle ang - this should only be used by the
c     dihedralEnergy function
c
c     (C) Christopher Woods
c
c##############################################################

      integer ipar,i,nterms
      double precision ang,ik1,ik2,ik3,ik4
      
      dihEnergy = ZERO
      
c     no invalid dihedral parameters!
      if (ipar.le.0 .or. ipar.gt.MAXDIHPARAM) return
      
      nterms = NDihTerms(ipar)
      if (nterms.le.0) return
      
      do i=1,nterms
        if (DihTermList(ipar,i).gt.0) then
          ik1 = DihTerm(DihTermList(ipar,i),1)
          ik2 = DihTerm(DihTermList(ipar,i),2)
          ik3 = DihTerm(DihTermList(ipar,i),3)
          ik4 = DihTerm(DihTermList(ipar,i),4)
          dihEnergy = dihEnergy + (ik1*(1.0 + ik2*(dcos(ik3*ang + ik4))))
        endif
      enddo
      
      return
      end
