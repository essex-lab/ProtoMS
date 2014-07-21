
      ipar1 = ZSolPar(funcisol1,funcat1)
      ipar2 = ZSolPar(funcisol2,funcat2)

      funcdx = dabs(ZSolCoords(funcisol1,funcat1,1)-ZSolCoords(funcisol2,funcat2,1))
      funcdy = dabs(ZSolCoords(funcisol1,funcat1,2)-ZSolCoords(funcisol2,funcat2,2))
      funcdz = dabs(ZSolCoords(funcisol1,funcat1,3)-ZSolCoords(funcisol2,funcat2,3))
      funcisvn1 = 0
      funcisvn2 = 0
      funciprot1 = 0
      funciprot2 = 0    
      funcdogeom = .false.
                     
      if (NSolVariable(funcisol1).gt.0 .or. NSolVariable(funcisol2).gt.0) then
        funcdogeom = .true.
        
        if (NSolVariable(funcisol1).gt.0) then
          funcitmp1 = NSolForwards(funcisol1)
        else
          funcitmp1 = funcisol1
        endif
        if (NSolVariable(funcisol2).gt.0) then
          funcitmp2 = NSolForwards(funcisol2)
        else
          funcitmp2 = funcisol2
        endif
        funcdxf = dabs(ZSolCoords(funcitmp1,funcat1,1)-ZSolCoords(funcitmp2,funcat2,1))
        funcdyf = dabs(ZSolCoords(funcitmp1,funcat1,2)-ZSolCoords(funcitmp2,funcat2,2))
        funcdzf = dabs(ZSolCoords(funcitmp1,funcat1,3)-ZSolCoords(funcitmp2,funcat2,3))
        
        if (NSolVariable(funcisol1).gt.0) then
          funcitmp1 = NSolBackwards(funcisol1)
        else
          funcitmp1 = funcisol1
        endif
        if (NSolVariable(funcisol2).gt.0) then
          funcitmp2 = NSolBackwards(funcisol2)
        else
          funcitmp2 = funcisol2
        endif
        funcdxb = dabs(ZSolCoords(funcitmp1,funcat1,1)-ZSolCoords(funcitmp2,funcat2,1))
        funcdyb = dabs(ZSolCoords(funcitmp1,funcat1,2)-ZSolCoords(funcitmp2,funcat2,2))
        funcdzb = dabs(ZSolCoords(funcitmp1,funcat1,3)-ZSolCoords(funcitmp2,funcat2,3))
      endif
