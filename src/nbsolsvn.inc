
      ipar1 = ZSolPar(funcisol1,funcat1)
      ipar2 = ZSvnPar(funcisvn2,funcat2)
      
      funcisol2 = 0
      funcisvn1 = 0
      funciprot1 = 0
      funciprot2 = 0       
      funcdogeom = .false.

      funcdx = dabs(ZSolCoords(funcisol1,funcat1,1)-ZSvnCoords(funcisvn2,funcat2,1))
      funcdy = dabs(ZSolCoords(funcisol1,funcat1,2)-ZSvnCoords(funcisvn2,funcat2,2))
      funcdz = dabs(ZSolCoords(funcisol1,funcat1,3)-ZSvnCoords(funcisvn2,funcat2,3))

      if (NSolVariable(funcisol1).gt.0) then
        funcdogeom = .true.  
      
        funcitmp = NSolForwards(funcisol1)
        funcdxf = dabs(ZSolCoords(funcitmp,funcat1,1)-ZSvnCoords(funcisvn2,funcat2,1))
        funcdyf = dabs(ZSolCoords(funcitmp,funcat1,2)-ZSvnCoords(funcisvn2,funcat2,2))
        funcdzf = dabs(ZSolCoords(funcitmp,funcat1,3)-ZSvnCoords(funcisvn2,funcat2,3))
        
        funcitmp = NSolBackwards(funcisol1)
        funcdxb = dabs(ZSolCoords(funcitmp,funcat1,1)-ZSvnCoords(funcisvn2,funcat2,1))
        funcdyb = dabs(ZSolCoords(funcitmp,funcat1,2)-ZSvnCoords(funcisvn2,funcat2,2))
        funcdzb = dabs(ZSolCoords(funcitmp,funcat1,3)-ZSvnCoords(funcisvn2,funcat2,3))
      endif
