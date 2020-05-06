
      ipar1 = ZPar(funciprot1,funcat1)
      ipar2 = ZSolPar(funcisol2,funcat2)

      funcdx = dabs(ZCoords(funciprot1,funcat1,1)-ZSolCoords(funcisol2,funcat2,1))
      funcdy = dabs(ZCoords(funciprot1,funcat1,2)-ZSolCoords(funcisol2,funcat2,2))
      funcdz = dabs(ZCoords(funciprot1,funcat1,3)-ZSolCoords(funcisol2,funcat2,3))
      funcisol1 = 0
      funciprot2 = 0
      funcisvn1 = 0
      funcisvn2 = 0
      funcdogeom = .false.
      
      if (NSolVariable(funcisol2).gt.0) then
        funcdogeom = .true.
        
        funcitmp1 = NSolForwards(funcisol2)
        funcdxf = dabs(ZCoords(funciprot1,funcat1,1)-ZSolCoords(funcitmp1,funcat2,1))
        funcdyf = dabs(ZCoords(funciprot1,funcat1,2)-ZSolCoords(funcitmp1,funcat2,2))
        funcdzf = dabs(ZCoords(funciprot1,funcat1,3)-ZSolCoords(funcitmp1,funcat2,3))
               
        funcitmp1 = NSolBackwards(funcisol2)
        funcdxb = dabs(ZCoords(funciprot1,funcat1,1)-ZSolCoords(funcitmp1,funcat2,1))
        funcdyb = dabs(ZCoords(funciprot1,funcat1,2)-ZSolCoords(funcitmp1,funcat2,2))
        funcdzb = dabs(ZCoords(funciprot1,funcat1,3)-ZSolCoords(funcitmp1,funcat2,3))
      endif
