
        ipar1 = ZPar(funciprot1,funcat1)
        ipar2 = ZPar(funciprot2,funcat2)

        funcdx = dabs(ZCoords(funciprot1,funcat1,1)-ZCoords(funciprot2,funcat2,1))
        funcdy = dabs(ZCoords(funciprot1,funcat1,2)-ZCoords(funciprot2,funcat2,2))
        funcdz = dabs(ZCoords(funciprot1,funcat1,3)-ZCoords(funciprot2,funcat2,3))
        funcisol1 = 0
        funcisol2 = 0
        funcisvn1 = 0
        funcisvn2 = 0
        funcdogeom = .false.
        
