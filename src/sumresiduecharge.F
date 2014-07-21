      double precision function sumResidueCharge(iprot,ires)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'proteins.inc'
      
      integer iprot,ires,i,ipar
      double precision ichg,charge
      
      sumResidueCharge = ZERO
      ichg = ZERO
      
c     sum up the backbone charge
      do i=1,4
        ipar = ZPar(iprot,BBZid(iprot,ires,i))
        ichg = ichg + charge(ipar)
      enddo
      
      do i=1,NSCAtoms(iprot,ires)
        ipar = ZPar(iprot,SCZid(iprot,ires,i))
        ichg = ichg + charge(ipar)
      enddo

      sumResidueCharge = ichg
      
      if (abs(sumResidueCharge).lt.0.00001) sumResidueCharge = ZERO
      
      return
      end
