
      ipar1 = ZGCSPar(funcigcs1,funcat1)
      ipar2 = ZGCSPar(funcigcs2,funcat2)
      funcdx = dabs(ZGCSCoords(funcigcs1,funcat1,1)-ZGCSCoords(funcigcs2,funcat2,1))
      funcdy = dabs(ZGCSCoords(funcigcs1,funcat1,2)-ZGCSCoords(funcigcs2,funcat2,2))
      funcdz = dabs(ZGCSCoords(funcigcs1,funcat1,3)-ZGCSCoords(funcigcs2,funcat2,3))
      funcisvn1 = 0
      funcisvn2 = 0
      funciprot1 = 0
      funciprot2 = 0    
      funcdogeom = .false.
      thet = thetaparam(funcigcs1)*thetaparam(funcigcs2)
