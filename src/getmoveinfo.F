      subroutine getMoveInfo()
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'splitter.inc'
      include 'constants.inc'
      include 'solvents.inc'
      include 'simparams.inc'
      include 'enums.inc'
      include 'moves.inc'
      include 'pressure.inc'
c########################################################################
c
c       This routine scans the words to find the optional data
c       that details move probabilities and how often to print
c       out move data
c
c       (C) Christopher Woods, October 2004
c            modified MSB Dec 08/09
c#########################################################################

      integer i,lgth,lineLength,getInt
       
      if (nwords.lt.3) return
                            
c     the format of these lines are;
c     ??? ??? (newprob) (protein=1) (solute=2) (solvent=3) (printinfo=10) #comment
      do i=3,nwords
        if (words(i)(1:1).eq.'#') goto 100
        lgth = lineLength(words(i))
        if (words(i).eq.'newprob') then
c         reset all of the probability data to 0 - nothing can move!
          pProtein = 0
          pSolvent = 0
          pSolute = 0  
          pTitrate = 0
          pTheta = 0
          pVolume = 0
          pRotamer = 0
          pInsertion = 0
          pDeletion = 0
          pGCSolute = 0
          pSample = 0
          pSolTheta = 0
        else if (words(i)(1:7).eq.'protein') then
          pProtein = getInt(words(i)(9:lgth),0)
        else if (words(i)(1:6).eq.'solute') then
          pSolute = getInt(words(i)(8:lgth),0)
        else if (words(i)(1:7).eq.'solvent') then
          pSolvent = getInt(words(i)(9:lgth),0)
        else if (words(i)(1:7).eq.'titrate') then
          pTitrate = getInt(words(i)(9:lgth),0)
        else if (words(i)(1:6).eq.'volume') then
          pVolume = getInt(words(i)(8:lgth),0)
        else if (words(i)(1:5).eq.'theta') then
          pTheta = getInt(words(i)(7:lgth),0)
        else if (words(i)(1:7).eq.'rotamer') then
          pRotamer = getInt(words(i)(9:lgth),0)
        else if (words(i)(1:9).eq.'insertion') then
          pInsertion = getInt(words(i)(11:lgth),0)
        else if (words(i)(1:8).eq.'deletion') then
          pDeletion = getInt(words(i)(10:lgth),0)
        else if (words(i)(1:8).eq.'gcsolute') then
          pGCSolute = getInt(words(i)(10:lgth),0)
        else if (words(i)(1:6).eq.'sample') then
          pSample = getInt(words(i)(8:lgth),0)
        else if (words(i)(1:8).eq.'soltheta') then
          pSolTheta = getInt(words(i)(10:lgth),0)
        else if (words(i)(1:9).eq.'printmove') then
          printMoveEvery = getInt(words(i)(11:lgth),0)
        
          if (printMoveEvery.lt.0) printMoveEvery = 0
          write(printstring,*) "Printing move info every ",printMoveEvery," steps"
          call printLine(INFO,printstring)
        endif
      enddo
100   continue      

 
c     eliminate invalid values
      if (NProteins.le.0) then 
        pProtein = 0
        pRotamer = 0
      endif

      if (pressure.eq.ZERO) pVolume = 0

      if (NSolvents.le.0) then
        pSolvent = 0
        pVolume = 0
        
       
      endif

      if (NSolutes.le.0) then
       pSolute = 0
      endif

      if (NGCSolutes.le.0) then
       pGCSolute = 0
      endif

      if (NThetaSolKind.le.0) pSolTheta = 0

      if (NProteins.le.0) pTitrate = 0

c     now calculate their values
      probTotal = pProtein+pSolute+pSolvent+pTitrate+pTheta+pVolume+pRotamer+pInsertion+pDeletion+pGCSolute+pSample+pSolTheta
      
c     if nothing is moving then set default values
      if (probTotal.eq.0) then
        call printLine(INFO,"Setting move probabilities to default values")
        pProtein = NProteins*20
        pSolute = NSolutes*20
        pSolvent = NSolvents*20
        pInsertion = 0
        pDeletion = 0
        pGCSolute = 0
        pTheta = 0
        pSample = 0
        pSolTheta = 0
        if (pressure .ne. ZERO) pVolume = 1
        probTotal = pProtein+pSolute+pSolvent+pTitrate+pTheta+pVolume+pRotamer+pInsertion+pDeletion+pGCSolute+pSolTheta
      endif
     
      probProtein = pProtein
      probSolute = pProtein+pSolute
      probSolvent = probSolute+pSolvent
      probTitrate = probSolvent+pTitrate
      probVolume = pVolume+probTitrate
      probTheta = pTheta+probVolume  
      probRotamer = pRotamer+probTheta
      probInsertion = pInsertion+probRotamer  
      probDeletion = pDeletion+probInsertion
      probGCSolute = pGCSolute+probDeletion
      probSample = pSample+probGCSolute
      probSolTheta = pSolTheta+probSample

      write(printstring,75) pProtein,pSolute,pSolvent,pInsertion
      call printLine(INFO,printstring)
      write(printstring,76) pTitrate,pTheta,pVolume,pDeletion
      call printLine(INFO,printstring)
      write(printstring,77) pRotamer,pSample,pGCsolute,pSolTheta
      call printLine(INFO,printstring)                   
                   
75    format('Move probabilies: Protein = ',i8,' Solute = ',i8,' Solvent = ',i8,' Insertion = ',i8)
76    format('                  Titrate = ',i8,' Theta  = ',i8,' Volume  = ',i8,' Deletion  = ',i8)
77    format('                  Rotamer = ',i8,' Sample = ',i8,' GCsolute = ',i8,' SolTheta = ',i8)

      end
