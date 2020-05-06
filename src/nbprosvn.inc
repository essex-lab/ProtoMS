
      ipar1 = ZPar(funciprot1,funcat1)
      ipar2 = ZSvnPar(funcisvn2,funcat2)

      funcdx = dabs(ZCoords(funciprot1,funcat1,1)-ZSvnCoords(funcisvn2,funcat2,1))
      funcdy = dabs(ZCoords(funciprot1,funcat1,2)-ZSvnCoords(funcisvn2,funcat2,2))
      funcdz = dabs(ZCoords(funciprot1,funcat1,3)-ZSvnCoords(funcisvn2,funcat2,3))
      funcisol1 = 0
      funcisol2 = 0
      funciprot2 = 0
      funcisvn1 = 0
      funcdogeom = .false.
