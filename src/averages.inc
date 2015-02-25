
c  This file has the variables used to store all of the energy averages

c       average energies (and average squares)
      double precision AvgProProC(MAXPROTEINS,MAXPROTEINS)
      double precision AvgProProC2(MAXPROTEINS,MAXPROTEINS)
      double precision AvgProProLJ(MAXPROTEINS,MAXPROTEINS)
      double precision AvgProProLJ2(MAXPROTEINS,MAXPROTEINS)
      double precision AvgProBnd(MAXPROTEINS),AvgProBnd2(MAXPROTEINS)
      double precision AvgProAng(MAXPROTEINS),AvgProAng2(MAXPROTEINS)
      double precision AvgProDih(MAXPROTEINS),AvgProDih2(MAXPROTEINS)
      double precision AvgProUBy(MAXPROTEINS),AvgProUBy2(MAXPROTEINS)
      double precision AvgProNBC(MAXPROTEINS),AvgProNBC2(MAXPROTEINS)
      double precision AvgProNBLJ(MAXPROTEINS),AvgProNBLJ2(MAXPROTEINS)
      double precision AvgSolSolC(MAXSOLUTES,MAXSOLUTES)
      double precision AvgSolSolC2(MAXSOLUTES,MAXSOLUTES)
      double precision AvgSolSolLJ(MAXSOLUTES,MAXSOLUTES)
      double precision AvgSolSolLJ2(MAXSOLUTES,MAXSOLUTES)
      double precision AvgSolBnd(MAXSOLUTES),AvgSolBnd2(MAXSOLUTES)
      double precision AvgSolAng(MAXSOLUTES),AvgSolAng2(MAXSOLUTES)
      double precision AvgSolDih(MAXSOLUTES),AvgSolDih2(MAXSOLUTES)
      double precision AvgSolUBy(MAXSOLUTES),AvgSolUBy2(MAXSOLUTES)
      double precision AvgSolNBC(MAXSOLUTES),AvgSolNBC2(MAXSOLUTES)
      double precision AvgSolNBLJ(MAXSOLUTES),AvgSolNBLJ2(MAXSOLUTES)
      double precision AvgProSolC(MAXPROTEINS,MAXSOLUTES)
      double precision AvgProSolC2(MAXPROTEINS,MAXSOLUTES)
      double precision AvgProSolLJ(MAXPROTEINS,MAXSOLUTES)
      double precision AvgProSolLJ2(MAXPROTEINS,MAXSOLUTES)
      double precision AvgSvnSvnC,AvgSvnSvnC2
      double precision AvgSvnSvnLJ,AvgSvnSvnLJ2
      double precision AvgSolSvnC(MAXSOLUTES),AvgSolSvnC2(MAXSOLUTES)
      double precision AvgSolSvnLJ(MAXSOLUTES),AvgSolSvnLJ2(MAXSOLUTES)
      double precision AvgProSvnC(MAXPROTEINS),AvgProSvnC2(MAXPROTEINS)
      double precision AvgProSvnLJ(MAXPROTEINS),AvgProSvnLJ2(MAXPROTEINS)
c       average forwards energies (and average squares)
      double precision AvgProProCF(MAXPROTEINS,MAXPROTEINS)
      double precision AvgProProCF2(MAXPROTEINS,MAXPROTEINS)
      double precision AvgProProLJF(MAXPROTEINS,MAXPROTEINS)
      double precision AvgProProLJF2(MAXPROTEINS,MAXPROTEINS)
      double precision AvgProBndF(MAXPROTEINS),AvgProBndF2(MAXPROTEINS)
      double precision AvgProAngF(MAXPROTEINS),AvgProAngF2(MAXPROTEINS)
      double precision AvgProDihF(MAXPROTEINS),AvgProDihF2(MAXPROTEINS)
      double precision AvgProUByF(MAXPROTEINS),AvgProUByF2(MAXPROTEINS)
      double precision AvgProNBCF(MAXPROTEINS),AvgProNBCF2(MAXPROTEINS)
      double precision AvgProNBLJF(MAXPROTEINS),AvgProNBLJF2(MAXPROTEINS)
      double precision AvgSolSolCF(MAXSOLUTES,MAXSOLUTES)
      double precision AvgSolSolCF2(MAXSOLUTES,MAXSOLUTES)
      double precision AvgSolSolLJF(MAXSOLUTES,MAXSOLUTES)
      double precision AvgSolSolLJF2(MAXSOLUTES,MAXSOLUTES)
      double precision AvgSolBndF(MAXSOLUTES),AvgSolBndF2(MAXSOLUTES)
      double precision AvgSolAngF(MAXSOLUTES),AvgSolAngF2(MAXSOLUTES)
      double precision AvgSolDihF(MAXSOLUTES),AvgSolDihF2(MAXSOLUTES)
      double precision AvgSolUByF(MAXSOLUTES),AvgSolUByF2(MAXSOLUTES)
      double precision AvgSolNBCF(MAXSOLUTES),AvgSolNBCF2(MAXSOLUTES)
      double precision AvgSolNBLJF(MAXSOLUTES),AvgSolNBLJF2(MAXSOLUTES)
      double precision AvgProSolCF(MAXPROTEINS,MAXSOLUTES)
      double precision AvgProSolCF2(MAXPROTEINS,MAXSOLUTES)
      double precision AvgProSolLJF(MAXPROTEINS,MAXSOLUTES)
      double precision AvgProSolLJF2(MAXPROTEINS,MAXSOLUTES)
      double precision AvgSvnSvnCF,AvgSvnSvnCF2
      double precision AvgSvnSvnLJF,AvgSvnSvnLJF2
      double precision AvgSolSvnCF(MAXSOLUTES),AvgSolSvnCF2(MAXSOLUTES)
      double precision AvgSolSvnLJF(MAXSOLUTES),AvgSolSvnLJF2(MAXSOLUTES)
      double precision AvgProSvnCF(MAXPROTEINS),AvgProSvnCF2(MAXPROTEINS)
      double precision AvgProSvnLJF(MAXPROTEINS),AvgProSvnLJF2(MAXPROTEINS)
c       average backwards energies (and average squares)
      double precision AvgProProCB(MAXPROTEINS,MAXPROTEINS)
      double precision AvgProProCB2(MAXPROTEINS,MAXPROTEINS)
      double precision AvgProProLJB(MAXPROTEINS,MAXPROTEINS)
      double precision AvgProProLJB2(MAXPROTEINS,MAXPROTEINS)
      double precision AvgProBndB(MAXPROTEINS),AvgProBndB2(MAXPROTEINS)
      double precision AvgProAngB(MAXPROTEINS),AvgProAngB2(MAXPROTEINS)
      double precision AvgProDihB(MAXPROTEINS),AvgProDihB2(MAXPROTEINS)
      double precision AvgProUByB(MAXPROTEINS),AvgProUByB2(MAXPROTEINS)
      double precision AvgProNBCB(MAXPROTEINS),AvgProNBCB2(MAXPROTEINS)
      double precision AvgProNBLJB(MAXPROTEINS),AvgProNBLJB2(MAXPROTEINS)
      double precision AvgSolSolCB(MAXSOLUTES,MAXSOLUTES)
      double precision AvgSolSolCB2(MAXSOLUTES,MAXSOLUTES)
      double precision AvgSolSolLJB(MAXSOLUTES,MAXSOLUTES)
      double precision AvgSolSolLJB2(MAXSOLUTES,MAXSOLUTES)
      double precision AvgSolBndB(MAXSOLUTES),AvgSolBndB2(MAXSOLUTES)
      double precision AvgSolAngB(MAXSOLUTES),AvgSolAngB2(MAXSOLUTES)
      double precision AvgSolDihB(MAXSOLUTES),AvgSolDihB2(MAXSOLUTES)
      double precision AvgSolUByB(MAXSOLUTES),AvgSolUByB2(MAXSOLUTES)
      double precision AvgSolNBCB(MAXSOLUTES),AvgSolNBCB2(MAXSOLUTES)
      double precision AvgSolNBLJB(MAXSOLUTES),AvgSolNBLJB2(MAXSOLUTES)
      double precision AvgProSolCB(MAXPROTEINS,MAXSOLUTES)
      double precision AvgProSolCB2(MAXPROTEINS,MAXSOLUTES)
      double precision AvgProSolLJB(MAXPROTEINS,MAXSOLUTES)
      double precision AvgProSolLJB2(MAXPROTEINS,MAXSOLUTES)
      double precision AvgSvnSvnCB,AvgSvnSvnCB2
      double precision AvgSvnSvnLJB,AvgSvnSvnLJB2
      double precision AvgSolSvnCB(MAXSOLUTES),AvgSolSvnCB2(MAXSOLUTES)
      double precision AvgSolSvnLJB(MAXSOLUTES),AvgSolSvnLJB2(MAXSOLUTES)
      double precision AvgProSvnCB(MAXPROTEINS),AvgProSvnCB2(MAXPROTEINS)
      double precision AvgProSvnLJB(MAXPROTEINS),AvgProSvnLJB2(MAXPROTEINS)
c     Bits for GC-solutes
      double precision AvgProGCSC(MAXPROTEINS),AvgProGCSC2(MAXPROTEINS)
      double precision AvgProGCSLJ(MAXPROTEINS),AvgProGCSLJ2(MAXPROTEINS)  
      double precision AvgProGCSCF(MAXPROTEINS),AvgProGCSCF2(MAXPROTEINS)
c      double precision AvgProGCSCLJF(MAXPROTEINS)
      double precision AvgProGCSLJF2(MAXPROTEINS),AvgProGCSLJF(MAXPROTEINS)
      double precision AvgProGCSCB(MAXPROTEINS),AvgProGCSCB2(MAXPROTEINS)
c      double precision AvgProGCSCLJB(MAXPROTEINS)
      double precision AvgProGCSLJB2(MAXPROTEINS),AvgProGCSLJB(MAXPROTEINS)
      
      double precision AvgSolGCSC(MAXSOLUTES),AvgSolGCSC2(MAXSOLUTES)
      double precision AvgSolGCSLJ(MAXSOLUTES),AvgSolGCSLJ2(MAXSOLUTES)  
      double precision AvgSolGCSCF(MAXSOLUTES),AvgSolGCSCF2(MAXSOLUTES)
      double precision AvgSolGCSLJF(MAXSOLUTES),AvgSolGCSLJB(MAXSOLUTES)      
c      double precision AvgSolGCSCLJF(MAXSOLUTES)
      double precision AvgSolGCSLJF2(MAXSOLUTES)
      double precision AvgSolGCSCB(MAXSOLUTES),AvgSolGCSCB2(MAXSOLUTES)
c      double precision AvgSolGCSCLJB(MAXSOLUTES)
      double precision AvgSolGCSLJB2(MAXSOLUTES)
      
      double precision AvgGCSSvnC(MAXSOLVENTS),AvgGCSSvnC2(MAXSOLVENTS)
      double precision AvgGCSSvnLJ(MAXSOLVENTS),AvgGCSSvnLJ2(MAXSOLVENTS)
      double precision AvgGCSSvnCF(MAXSOLVENTS),AvgGCSSvnCF2(MAXSOLVENTS)
c      double precision AvgGCSSvnCLJF(MAXSOLVENTS)
      double precision AvgGCSSvnLJF(MAXSOLVENTS),AvgGCSSvnLJB(MAXSOLVENTS)
      double precision AvgGCSSvnLJF2(MAXSOLVENTS),AvgGCSSvnLJB2(MAXSOLVENTS)
      double precision AvgGCSSvnCB(MAXSOLVENTS),AvgGCSSvnCB2(MAXSOLVENTS)
c      double precision AvgGCSSvnCLJB(MAXSOLVENTS)
      double precision AvgGCSGCSC,AvgGCSGCSC2,AvgIndGCSGCSCF
      double precision AvgGCSGCSLJ,AvgGCSGCSLJ2
      double precision AvgGCSGCSCF,AvgGCSGCSCF2,AvgGCSC2
      double precision AvgGCSGCSLJF2,AvgIndGCSGCSC,AvgGCSCB,AvgGCSGCSLJF
      double precision AvgGCSGCSCB,AvgGCSGCSCB2,AvgGCSLJ,AvgGCSCF2
      double precision AvgGCSGCSLJB
      double precision AvgGCSGCSLJB2,AvgGCSCB2,AvgGCSLJF2,AvgGCSLJB2, AvgGCSLJ2
      
      
      
      common /DAVGENERGY1/ AvgProProC,AvgProProC2,AvgProProLJ,AvgProProLJ2,
     .     AvgProBnd,AvgProBnd2,AvgProAng,AvgProAng2,AvgProDih,AvgProDih2,
     .     AvgProUby,AvgProUby2,AvgSolUBy,AvgSolUBy2,
     .     AvgProNBC,AvgProNBC2,AvgProNBLJ,AvgProNBLJ2,AvgSolSolC,
     .     AvgSolSolC2,AvgSolSolLJ,AvgSolSolLJ2,AvgSolBnd,AvgSolBnd2,
     .     AvgSolAng,AvgSolAng2,AvgSolDih,AvgSolDih2,AvgSolNBC,AvgSolNBC2,
     .     AvgSolNBLJ,AvgSolNBLJ2,AvgProSolC,AvgProSolC2,AvgProSolLJ,AvgProSolLJ2,
     .     AvgSvnSvnC,AvgSvnSvnC2,AvgSvnSvnLJ,AvgSvnSvnLJ2,AvgSolSvnC,AvgSolSvnC2,
     .     AvgSolSvnLJ,AvgSolSvnLJ2,AvgProSvnC,AvgProSvnC2,AvgProSvnLJ,AvgProSvnLJ2
      
      common /DAVGENERGY2/ avgtotalenergyb2,avgtotalenergyb,avgtotalenergyf2,avgtotalenergyf,
     .     avgprosvnljb2,avgprosvnljb,avgprosvncb2,avgprosvncb,
     .     avgsolsvnljb2,avgsolsvnljb,avgsolsvncb2,avgsolsvncb,avgsvnsvnljb2,
     .     avgsvnsvnljb,avgsvnsvncb2,avgsvnsvncb,avgprosolljb2,avgprosolljb,
     .     avgprosolcb2,avgprosolcb,avgsolnbljb2,avgsolnbljb,avgsolnbcb2,
     .     avgsolnbcb,avgsoldihb2,avgsoldihb,avgsolangb2,avgsolangb,
     .     avgsolbndb2,avgsolbndb,avgsolsolljb2,avgsolsolljb,avgsolsolcb2,
     .     avgsolsolcb,avgpronbljb2,avgpronbljb,avgpronbcb2,avgpronbcb,
     .     avgprodihb2,avgprodihb,avgproangb2,avgproangb,avgprobndb2,
     .     avgproubyb,avgproubyb2,avgproubyf,avgproubyf2,
     .     avgprobndb,avgproproljb2,avgproproljb,avgproprocb2,avgproprocb,
     .     avgprosvnljf2,avgprosvnljf,avgprosvncf2,avgprosvncf,avgsolsvnljf2,
     .     avgsolsvnljf,avgsolsvncf2,avgsolsvncf,avgsvnsvnljf2,avgsvnsvnljf,
     .     avgsvnsvncf2,avgsvnsvncf,avgprosolljf2,avgprosolljf,avgprosolcf2,
     .     avgprosolcf,avgsolnbljf2,avgsolnbljf,avgsolnbcf2,avgsolnbcf,
     .     avgsoldihf2,avgsoldihf,avgsolangf2,avgsolangf,avgsolbndf2,
     .     avgsolubyf,avgsolubyf2,avgsolubyb,avgsolubyb2,
     .     avgsolbndf,avgsolsolljf2,avgsolsolljf,avgsolsolcf2,avgsolsolcf,
     .     avgpronbljf2,avgpronbljf,avgpronbcf2,avgpronbcf,avgprodihf2,
     .     avgprodihf,avgproangf2,avgproangf,avgprobndf2,avgprobndf,
     .     avgproproljf2,avgproproljf,avgproprocf2,avgproprocf
     
      common /DAVGENERGY3/ AvgProGCSC,AvgProGCSCF,AvgProGCSCB,AvgProGCSLJF2,AvgProGCSLJB2,
     .     AvgProGCSLJ,AvgProGCSC2,AvgProGCSCF2,AvgProGCSCB2,AvgProGCSLJF,AvgProGCSLJB,
     .     AvgProGCSLJ2,AvgGCSSvnC,AvgGCSSvnCF,AvgGCSSvnCB,AvgGCSGCSLJF2,AvgGCSGCSLJF,AvgGCSGCSLJB2,AvgGCSGCSLJB,
     .     AvgGCSSvnLJ,AvgGCSSvnC2,AvgGCSSvnCF2,AvgGCSSvnCB2,AvgGCSSvnLJF2,AvgGCSSvnLJB2,AvgGCSSvnLJF,
     .     AvgGCSSvnLJ2,AvgIndGCSGCSC,AvgIndGCSGCSCF,AvgGCSCB,AvgGCSGCSLJ2,AvgGCSGCSLJ,AvgGCSSvnLJB,
     .     AvgGCSLJ,AvgGCSC2,AvgGCSCF2,AvgGCSCB2,AvgGCSLJF2,AvgGCSLJB2,AvgGCSGCSC,AvgGCSGCSC2,
     .     AvgGCSLJ2, AvgSolGCSC,AvgSolGCSCF,AvgSolGCSCB,AvgGCSGCSCF2,AvgGCSGCSCF,AvgGCSGCSCB2,AvgGCSGCSCB,
     .     AvgSolGCSLJ,AvgSolGCSC2,AvgSolGCSCF2,AvgSolGCSCB2,AvgSolGCSLJF2,AvgSolGCSLJB2,
     .     AvgSolGCSLJ2, AvgSolGCSLJF, AvgSolGCSLJB

      

       save /DAVGENERGY1/
       save /DAVGENERGY2/
       save /DAVGENERGY3/
       
                            
c     average energies / free energies
      double precision AvgTotalEnergy,AvgTotalEnergy2
      double precision AvgTotalEnergyF,AvgTotalEnergyF2
      double precision AvgTotalEnergyB,AvgTotalEnergyB2
      double precision AvgBackEnergy,AvgForwardsEnergy
      double precision AvgBackEnergy1,AvgForwardsEnergy1
      double precision AvgBackEnergy2,AvgForwardsEnergy2
c     SG PME energies
      double precision AvgPmeEnergy,AvgPmeEnergyF,AvgPmeEnergyB
      integer NSteps
                    
      common /FREEENERGIES/ AvgTotalEnergy,AvgBackEnergy,AvgForwardsEnergy,
     .            AvgTotalEnergy2,AvgBackEnergy2,AvgForwardsEnergy2,
     .            AvgBackEnergy1,AvgForwardsEnergy1,AvgPmeEnergy,AvgPmeEnergyF,AvgPmeEnergyB,NSteps

      save /FREEENERGIES/
      
