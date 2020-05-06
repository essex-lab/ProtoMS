
      ipar1 = ZGCSPar(funcigcs1,funcat1)
      ipar2 = ZSvnPar(funcisvn2,funcat2)
      
      funcigcs2 = 0
      funcisvn1 = 0
      funciprot1 = 0
      funciprot2 = 0       
      funcdogeom = .false.
      thet = thetaparam(funcigcs1)
      funcdx = dabs(ZGCSCoords(funcigcs1,funcat1,1)-ZSvnCoords(funcisvn2,funcat2,1))
      funcdy = dabs(ZGCSCoords(funcigcs1,funcat1,2)-ZSvnCoords(funcisvn2,funcat2,2))
      funcdz = dabs(ZGCSCoords(funcigcs1,funcat1,3)-ZSvnCoords(funcisvn2,funcat2,3))

c      if (NGCSVariable(funcigcs1).gt.0) then
c        funcdogeom = .true.  
c      
c        funcitmp = NGCSForwards(funcigcs1)
c        funcdxf = dabs(ZGCSCoords(funcitmp,funcat1,1)-ZSvnCoords(funcisvn2,funcat2,1))
c        funcdyf = dabs(ZGCSCoords(funcitmp,funcat1,2)-ZSvnCoords(funcisvn2,funcat2,2))
c        funcdzf = dabs(ZGCSCoords(funcitmp,funcat1,3)-ZSvnCoords(funcisvn2,funcat2,3))
c        
c        funcitmp = NGCSBackwards(funcigcs1)
c        funcdxb = dabs(ZGCSCoords(funcitmp,funcat1,1)-ZSvnCoords(funcisvn2,funcat2,1))
c        funcdyb = dabs(ZGCSCoords(funcitmp,funcat1,2)-ZSvnCoords(funcisvn2,funcat2,2))
c        funcdzb = dabs(ZGCSCoords(funcitmp,funcat1,3)-ZSvnCoords(funcisvn2,funcat2,3))
c      endif
