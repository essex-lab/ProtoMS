
      subroutine printAverages
      include 'dimensions.inc'
      include 'metropolis.inc'
      include 'printer.inc'
      include 'flags.inc'
      include 'stacks.inc'
      include 'proteins.inc'
      include 'simparams.inc'
      include 'parameters.inc'
      include 'solvents.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'averages.inc'
      include 'enums.inc'
      include 'pressure.inc'
      include 'moves.inc'
      include 'energies.inc'
      include 'parallel.inc'
c###################################################################
c
c     This routine prints the averages, as of the current time!
c
c     (C) Christopher Woods          9 August 2002
c
c####################################################################

       double precision gfor,gback,stdDev
       double precision std,stdf,stdb,totalEnergy,nrgb,nrgf
       integer lineLength,lgth,lgth2,igcs,savenrg,Molcount
       integer i,j,isvn,itry,iacp,iprcnt,ires,X,Y,Z,s,t,initial


       if (DifferentLambdas) then
         write(printstring,10) lambdab,lambda,lambdaf
10       format("RESULTS for lambda ",f8.5," <-- ",f8.5," --> ",f8.5)
       else 
         write(printstring,11) lambda
11       format("RESULTS for lambda = ",f8.5)
       endif
       call prettyLine(RESULTS,printstring)
            
c       write(printstring,*) "Total Energy = ",totalEnergy(savenrg,nrgb,nrgf)
c       call printLine(RESULTS,printstring)


       write(printstring,"('Number of data steps = ',i8)") NSTEPS 
       call printLine(RESULTS,printstring)

       if (doLambdaRE) then
         write(printstring,"('Lambda replica = ',i8)") mylamid+1
         call printLine(RESULTS,printstring)
       endif

       if (doTemperatureRE) then
         write(printstring,"('Temperature replica = ',i8)") mytemid+1
         call printLine(RESULTS,printstring)
       endif

       if (doMultiGCMC) then
         write(printstring,"('GCMC replica = ',i8)") mygcid+1
         call printLine(RESULTS,printstring)
       endif

       if (doTemperatureRE.and.doLambdaRE .or. doMultiGCMC.and.doLambdaRE) then
         write(printstring,"('Global replica = ',i8)") myid+1
         call printLine(RESULTS,printstring)
       endif

       write(printstring,"('Temperature = ',f12.4,' C')") temperatureC
       call printLine(RESULTS,printstring)

       if (doREST) then
         write(printstring,"('Effective temperature = ',f12.4,' C')") resttempC
         call printLine(RESULTS,printstring)
       endif

       write(printstring,"('Solvents,Proteins,GC-solutes = ',3i8)") NSolvents ,NProteins ,NGCSolutes
       call printLine(RESULTS,printstring)

       write(printstring,"('Simulation B value = ',f12.4)") B
       call printLine(RESULTS,printstring)

       write(printstring,"('Molecules in grid  = ',i8)") GCSCavity
       call printLine(RESULTS,printstring)

        if (SolventBox) then
         write(printstring,"('Pressure = ',f12.4,' atm')") pressure
         call printLine(RESULTS,printstring)
       
         write(printstring,"('Volume = ',f12.4,' A^3')") volume
         call printLine(RESULTS,printstring)
       
         std = stdDev(AvgVolume2,AvgVolume)
         write(printstring,"('Average volume = ',f12.4,' (',f12.4,')')") AvgVolume,std
         call printLine(RESULTS,printstring)
       endif

c      write out the random number seed used to generate this trajectory
       write(printstring,"('Random number seed = ',i12)") ranseed
       call printLine(RESULTS,printstring)

       if (NSteps.le.0) then
         write(printstring,*) "No data collected, thus no averages!"
         call printLine(RESULTS,printstring)
         return
       endif
                                  
       if (DifferentLambdas) then
         gfor = -dlog(AVGFORWARDSENERGY) / mt_beta_T
         gback = -dlog(AVGBACKENERGY) / mt_beta_T

         if (abs(gfor).lt.TINY) gfor = ZERO
         if (abs(gback).lt.TINY) gback = ZERO
       
         std = stdDev(AVGBACKENERGY2,AVGBACKENERGY1)
         write(printstring,19) "Backwards",gback,std
         call printLine(RESULTS,printstring)
         
         std = stdDev(AVGFORWARDSENERGY2,AVGFORWARDSENERGY1)
         write(printstring,19) "Forwards ",gfor,std
         call printLine(RESULTS,printstring)

19       format(a9,' Free Energy = ',f15.4,' ( ',f8.3,' )')
       endif
       write(printstring,*) "Averages. Stddev in brackets. Average | forwards backwards"
       call prettyLine(RESULTS,printstring)
       
       std = stdDev(AVGTOTALENERGY2,AVGTOTALENERGY)
       stdf = stdDev(AVGTOTALENERGYF2,AVGTOTALENERGYF)
       stdb = stdDev(AVGTOTALENERGYB2,AVGTOTALENERGYB)
       write(printstring,43) AVGTOTALENERGY,std,AVGTOTALENERGYF,stdf,AVGTOTALENERGYB,stdb
43     format('Total Energy ',f15.4,' ( ',f8.3,' ) | ',f15.4,' ( ',f8.3,' ) ',f15.4,' ( ',f8.3,' )')
       call printLine(RESULTS,printstring)

       if (NProteins.gt.0) then
         do i=1,NProteins
           if (NBonds(i).gt.0) then
c           internal bond energy
             std = stdDev(AvgProBnd2(i),AvgProBnd(i))
c             PRINT*, "std" ,std
             stdf = stdDev(AvgProBndF2(i),AvgProBndF(i))
             stdb = stdDev(AvgProBndB2(i),AvgProBndB(i))
             lgth = lineLength(ProteinName(i))
             write(printstring,"('Internal energy averages for protein ',i1,1x,a)") i,ProteinName(i)(1:lgth)
             call prettyLine(RESULTS,printstring)
             write(printstring,20) "BND",AvgProBnd(i),std,AvgProBndF(i),stdf,AvgProBndB(i),stdb
20           format(a3,1x,f15.4,' ( ',f8.3,' ) | ',f15.4,' ( ',f8.3,' ) ',f15.4,' ( ',f8.3,' )')
             call printLine(RESULTS,printstring)
             write(printstring,21) AvgProBndF(i)-AvgProBnd(i),AvgProBndB(i)-AvgProBnd(i)
             call printLine(RESULTS,printstring)
21           format('DIFF',31x,f15.4,14x,f15.4)
           endif
           
           if (NAngles(i).gt.0) then
c            internal angle energy
             std = stdDev(AvgProAng2(i),AvgProAng(i))
             stdf = stdDev(AvgProAngF2(i),AvgProAngF(i))
             stdb = stdDev(AvgProAngB2(i),AvgProAngB(i))
             write(printstring,20) "ANG",AvgProAng(i),std,AvgProAngF(i),stdf,AvgProAngB(i),stdb
             call printLine(RESULTS,printstring)
             write(printstring,21) AvgProAngF(i)-AvgProAng(i),AvgProAngB(i)-AvgProAng(i)
             call printLine(RESULTS,printstring)
           endif
                
           if (NDihedrals(i).gt.0) then               
c            internal dihedral energy
             std = stdDev(AvgProDih2(i),AvgProDih(i))
             stdf = stdDev(AvgProDihF2(i),AvgProDihF(i))
             stdb = stdDev(AvgProDihB2(i),AvgProDihB(i))
             write(printstring,20) "DIH",AvgProDih(i),std,AvgProDihF(i),stdf,AvgProDihB(i),stdb
             call printLine(RESULTS,printstring)
             write(printstring,21) AvgProDihF(i)-AvgProDih(i),AvgProDihB(i)-AvgProDih(i)
             call printLine(RESULTS,printstring)
           endif

           if (NUreyBradleys(i).gt.0) then               
c            internal Urey-Bradley energy
             std = stdDev(AvgProUBy2(i),AvgProUBy(i))
             stdf = stdDev(AvgProUByF2(i),AvgProUByF(i))
             stdb = stdDev(AvgProUByB2(i),AvgProUByB(i))
             write(printstring,20) "UBY",AvgProUBy(i),std,AvgProUByF(i),stdf,AvgProUByB(i),stdb
             call printLine(RESULTS,printstring)
             write(printstring,21) AvgProUByF(i)-AvgProUBy(i),AvgProUByB(i)-AvgProUBy(i)
             call printLine(RESULTS,printstring)
           endif

c          internal nbc energy
           std = stdDev(AvgProNBC2(i),AvgProNBC(i))
           stdf = stdDev(AvgProNBCF2(i),AvgProNBCF(i))
           stdb = stdDev(AvgProNBCB2(i),AvgProNBCB(i))
           write(printstring,20) "NBC",AvgProNBC(i),std,AvgProNBCF(i),stdf,AvgProNBCB(i),stdb
           call printLine(RESULTS,printstring)
           write(printstring,21) AvgProNBCF(i)-AvgProNBC(i),AvgProNBCB(i)-AvgProNBC(i)
           call printLine(RESULTS,printstring)

c          internal nblj energy
           std = stdDev(AvgProNBLJ2(i),AvgProNBLJ(i))
           stdf = stdDev(AvgProNBLJF2(i),AvgProNBLJF(i))
           stdb = stdDev(AvgProNBLJB2(i),AvgProNBLJB(i))
           write(printstring,20) "NBL",AvgProNBLJ(i),std,AvgProNBLJF(i),stdf,AvgProNBLJB(i),stdb
           call printLine(RESULTS,printstring)
           write(printstring,21) AvgProNBLJF(i)-AvgProNBLJ(i),AvgProNBLJB(i)-AvgProNBLJ(i)
           call printLine(RESULTS,printstring)
         enddo
       else
c       PRINT*, "Hello"
       endif           

       if (NSolutes.gt.0) then
         do i=1,NSolutes
c          skip duplicate solutes used for geometry variations / multicopy sampling
           if (.not.IgnoreSolute(i)) then
              
             lgth = lineLength(SoluteName(i))
             write(printstring,23) i,SoluteName(i)(1:lgth)
23           format('Internal energy averages for solute ',i4,' ',a)
             call prettyLine(RESULTS,printstring)
           
             if (NSolBonds(i).gt.0) then
c              bond
               std = stdDev(AvgSolBnd2(i),AvgSolBnd(i))
               stdf = stdDev(AvgSolBndf2(i),AvgSolBndf(i))
               stdb = stdDev(AvgSolBndb2(i),AvgSolBndb(i))
               write(printstring,20) 'BND',AvgSolBnd(i),std,AvgSolBndF(i),stdf,AvgSolBndB(i),stdb
               call printLine(RESULTS,printstring)
               write(printstring,21) AvgSolBndF(i)-AvgSolBnd(i),AvgSolBndB(i)-AvgSolBnd(i)
               call printLine(RESULTS,printstring)
             endif
           
             if (NSolAngles(i).gt.0) then
c              angle
               std = stdDev(AvgSolAng2(i),AvgSolAng(i))
               stdf = stdDev(AvgSolAngf2(i),AvgSolAngf(i))
               stdb = stdDev(AvgSolAngb2(i),AvgSolAngb(i))
               write(printstring,20) 'ANG',AvgSolAng(i),std,AvgSolAngF(i),stdf,AvgSolAngB(i),stdb
               call printLine(RESULTS,printstring)
               write(printstring,21) AvgSolAngF(i)-AvgSolAng(i),AvgSolAngB(i)-AvgSolAng(i)
               call printLine(RESULTS,printstring)
             endif

             if (NSolDihedrals(i).gt.0) then
c              dihedral
               std = stdDev(AvgSolDih2(i),AvgSolDih(i))
               stdf = stdDev(AvgSolDihf2(i),AvgSolDihf(i))
               stdb = stdDev(AvgSolDihb2(i),AvgSolDihb(i))
               write(printstring,20) 'DIH',AvgSolDih(i),std,AvgSolDihF(i),stdf,AvgSolDihB(i),stdb
               call printLine(RESULTS,printstring)
               write(printstring,21) AvgSolDihF(i)-AvgSolDih(i),AvgSolDihB(i)-AvgSolDih(i)
               call printLine(RESULTS,printstring)
             endif

             if (NSolUreyBradleys(i).gt.0) then
c              Urey-Bradley
               std = stdDev(AvgSolUBy2(i),AvgSolUBy(i))
               stdf = stdDev(AvgSolUByf2(i),AvgSolUByf(i))
               stdb = stdDev(AvgSolUByb2(i),AvgSolUByb(i))
               write(printstring,20) 'UBY',AvgSolUBy(i),std,AvgSolUByF(i),stdf,AvgSolUByB(i),stdb
               call printLine(RESULTS,printstring)
               write(printstring,21) AvgSolUByF(i)-AvgSolUBy(i),AvgSolUByB(i)-AvgSolUBy(i)
               call printLine(RESULTS,printstring)
             endif

c            nbc
             std = stdDev(AvgSolNBC2(i),AvgSolNBC(i))
             stdf = stdDev(AvgSolNBCf2(i),AvgSolNBCf(i))
             stdb = stdDev(AvgSolNBCb2(i),AvgSolNBCb(i))
             write(printstring,20) 'NBC',AvgSolNBC(i),std,AvgSolNBCF(i),stdf,AvgSolNBCB(i),stdb
             call printLine(RESULTS,printstring)
             write(printstring,21) AvgSolNBCF(i)-AvgSolNBC(i),AvgSolNBCB(i)-AvgSolNBC(i)
             call printLine(RESULTS,printstring)

c            nblj
             std = stdDev(AvgSolNBLJ2(i),AvgSolNBLJ(i))
             stdf = stdDev(AvgSolNBLJf2(i),AvgSolNBLJf(i))
             stdb = stdDev(AvgSolNBLJb2(i),AvgSolNBLJb(i))
             write(printstring,20) 'NBL',AvgSolNBLJ(i),std,AvgSolNBLJF(i),stdf,AvgSolNBLJB(i),stdb
             call printLine(RESULTS,printstring)
             write(printstring,21) AvgSolNBLJF(i)-AvgSolNBLJ(i),AvgSolNBLJB(i)-AvgSolNBLJ(i)
             call printLine(RESULTS,printstring)
           endif
         enddo
       endif
       
       if (NSolvents.gt.0) then
         call prettyLine(RESULTS,"Average solvent-solvent energy")
c        coulomb energy         
c         PRINT*, "av" ,AvgSvnSvnC2
         std = stdDev(AvgSvnSvnC2,AvgSvnSvnC)
         stdf = stdDev(AvgSvnSvnCF2,AvgSvnSvnCF)
         stdb = stdDev(AvgSvnSvnCB2,AvgSvnSvnCB)
         write(printstring,20) 'COU',AvgSvnSvnC,std,AvgSvnSvnCF,stdf,AvgSvnSvnCB,stdb
         call printLine(RESULTS,printstring)
         write(printstring,21) AvgSvnSvnCF-AvgSvnSvnC,AvgSvnSvnCB-AvgSvnSvnC
         call printLine(RESULTS,printstring)
         
c        lj energy
         std = stdDev(AvgSvnSvnLJ2,AvgSvnSvnLJ)
         stdf = stdDev(AvgSvnSvnLJF2,AvgSvnSvnLJF)
         stdb = stdDev(AvgSvnSvnLJB2,AvgSvnSvnLJB)
         write(printstring,20) 'LJ ',AvgSvnSvnLJ,std,AvgSvnSvnLJF,stdf,AvgSvnSvnLJB,stdb
         call printLine(RESULTS,printstring)
         write(printstring,21) AvgSvnSvnLJF-AvgSvnSvnLJ,AvgSvnSvnLJB-AvgSvnSvnLJ
         call printLine(RESULTS,printstring)
       endif
       
       if (NGCSolutes.gt.0) then
         call prettyLine(RESULTS,"Average GCS-GCS energy")
c        coulomb energy         
c         PRINT*, "av" ,AvgGCSGCSC2
         std = stdDev(AvgGCSGCSC2,AvgGCSGCSC)
         stdf = stdDev(AvgGCSGCSCF2,AvgGCSGCSCF)
         stdb = stdDev(AvgGCSGCSCB2,AvgGCSGCSCB)
         write(printstring,20) 'COU',AvgGCSGCSC,std,AvgGCSGCSCF,stdf,AvgGCSGCSCB,stdb
         call printLine(RESULTS,printstring)
         write(printstring,21) AvgGCSGCSCF-AvgGCSGCSC,AvgGCSGCSCB-AvgGCSGCSC
         call printLine(RESULTS,printstring)
         
c        lj energy
         std = stdDev(AvgGCSGCSLJ2,AvgGCSGCSLJ)
         stdf = stdDev(AvgGCSGCSLJF2,AvgGCSGCSLJF)
         stdb = stdDev(AvgGCSGCSLJB2,AvgGCSGCSLJB)
         write(printstring,20) 'LJ ',AvgGCSGCSLJ,std,AvgGCSGCSLJF,stdf,AvgGCSGCSLJB,stdb
         call printLine(RESULTS,printstring)
         write(printstring,21) AvgGCSGCSLJF-AvgGCSGCSLJ,AvgGCSGCSLJB-AvgGCSGCSLJ
         call printLine(RESULTS,printstring)
       endif
                     
         write(printstring,*) "Average inter-molecular interaction energies"
         call prettyLine(RESULTS,printstring)


       if (NProteins.gt.0 .and. NGCSolutes.gt.0) then
         do i=1,NProteins
           lgth = lineLength(ProteinName(i))
           write(printstring,933) i,ProteinName(i)(1:lgth)
933        format('Average protein-GCS energy for ',i4,' ',a)      
           call prettyLine(RESULTS,printstring)

c            coulomb
c             PRINT*, "en" ,AvgProGCSC
             std = stdDev(AvgProGCSC2(i),AvgProGCSC(i))
             stdf = stdDev(AvgProGCSCF2(i),AvgProGCSCF(i))
             stdb = stdDev(AvgProGCSCB2(i),AvgProGCSCB(i))
             write(printstring,20) 'COU',AvgProGCSC(i),std,AvgProGCSCF(i),stdf,
     .                                                      AvgProGCSCB(i),stdb
             call printLine(RESULTS,printstring)
             write(printstring,21) AvgProGCSCF(i)-AvgProGCSC(i),AvgProGCSCB(i)-AvgProGCSC(i)
             call printLine(RESULTS,printstring)

c            lj
             std = stdDev(AvgProGCSLJ2(i),AvgProGCSLJ(i))
             stdf = stdDev(AvgProGCSLJF2(i),AvgProGCSLJF(i))
             stdb = stdDev(AvgProGCSLJB2(i),AvgProGCSLJB(i))
             write(printstring,20) 'LJ',AvgProGCSLJ(i),std,AvgProGCSLJF(i),stdf,
     .                                                    AvgProGCSLJB(i),stdb
             call printLine(RESULTS,printstring)
             write(printstring,21) AvgProGCSLJF(i)-AvgProGCSLJ(i),AvgProGCSLJB(i)-AvgProGCSLJ(i)
             call printLine(RESULTS,printstring)
         enddo
       endif

       if (NGCSolutes.gt.0 .and. NSolvents.gt.0) then
         do i=1,NGCSolutes
           lgth = lineLength(GCSName(i))
           write(printstring,932) i,GCSName(i)(1:lgth)
932        format('Average GCS-solvent energy for ',i4,' ',a)      
           call prettyLine(RESULTS,printstring)

c            coulomb
c             PRINT*, "en" ,AvgGCSSvnC
             std = stdDev(AvgGCSSvnC2(i),AvgGCSSvnC(i))
             stdf = stdDev(AvgGCSSvnCF2(i),AvgGCSSvnCF(i))
             stdb = stdDev(AvgGCSSvnCB2(i),AvgGCSSvnCB(i))
             write(printstring,20) 'COU',AvgGCSSvnC(i),std,AvgGCSSvnCF(i),stdf,
     .                                                      AvgGCSSvnCB(i),stdb
             call printLine(RESULTS,printstring)
             write(printstring,21) AvgGCSSvnCF(i)-AvgGCSSvnC(i),AvgGCSSvnCB(i)-AvgGCSSvnC(i)
             call printLine(RESULTS,printstring)

c            lj
             std = stdDev(AvgGCSSvnLJ2(i),AvgGCSSvnLJ(i))
             stdf = stdDev(AvgGCSSvnLJF2(i),AvgGCSSvnLJF(i))
             stdb = stdDev(AvgGCSSvnLJB2(i),AvgGCSSvnLJB(i))
             write(printstring,20) 'LJ',AvgGCSSvnLJ(i),std,AvgGCSSvnLJF(i),stdf,
     .                                                    AvgGCSSvnLJB(i),stdb
             call printLine(RESULTS,printstring)
             write(printstring,21) AvgGCSSvnLJF(i)-AvgGCSSvnLJ(i),AvgGCSSvnLJB(i)-AvgGCSSvnLJ(i)
             call printLine(RESULTS,printstring)
         enddo
       endif

       
          
c"

        

         do i=1,NProteins-1
           do j=i+1,NProteins
             write(printstring,30) i,ProteinName(i)(1:20),j,ProteinName(j)(1:20)
30           format('Average protein-protein energy for ',i4,' ',a20,' - ',i4,' ',a20)             
             call prettyLine(RESULTS,printstring)

c            coulomb
             std = stdDev(AvgProProC2(i,j),AvgProProC(i,j))
             stdf = stdDev(AvgProProCF2(i,j),AvgProProCF(i,j))
             stdb = stdDev(AvgProProCB2(i,j),AvgProProCB(i,j))
             write(printstring,20) 'COU',AvgProProC(i,j),std,AvgProProCF(i,j),stdf,
     .                                                      AvgProProCB(i,j),stdb
             call printLine(RESULTS,printstring)
             write(printstring,21) AvgProProCF(i,j)-AvgProProC(i,j),AvgProProCB(i,j)-AvgProProC(i,j)
             call printLine(RESULTS,printstring)

c            lj
             std = stdDev(AvgProProLJ2(i,j),AvgProProLJ(i,j))
             stdf = stdDev(AvgProProLJF2(i,j),AvgProProLJF(i,j))
             stdb = stdDev(AvgProProLJB2(i,j),AvgProProLJB(i,j))
             write(printstring,20) 'LJ ',AvgProProLJ(i,j),std,AvgProProLJF(i,j),stdf,
     .                                                      AvgProProLJB(i,j),stdb
             call printLine(RESULTS,printstring)
             write(printstring,21) AvgProProLJF(i,j)-AvgProProLJ(i,j),AvgProProLJB(i,j)-AvgProProLJ(i,j)
             call printLine(RESULTS,printstring)
             
           enddo
         enddo
c        endif
       
       if (NProteins.gt.0 .and. NSolutes.gt.0) then
         do i=1,NProteins
           do j=1,NSolutes
             if (.not.IgnoreSolute(j)) then
               lgth = lineLength(ProteinName(i))
               lgth2 = lineLength(SoluteName(j))
               write(printstring,32) i,ProteinName(i)(1:lgth),j,SoluteName(j)(1:lgth2)
32             format('Average solute-protein energy for ',i4,' ',a,' - ',i4,' ',a)             
               call prettyLine(RESULTS,printstring)

c              coulomb
               std = stdDev(AvgProSolC2(i,j),AvgProSolC(i,j))
               stdf = stdDev(AvgProSolCF2(i,j),AvgProSolCF(i,j))
               stdb = stdDev(AvgProSolCB2(i,j),AvgProSolCB(i,j))
               write(printstring,20) 'COU',AvgProSolC(i,j),std,AvgProSolCF(i,j),stdf,
     .                                                      AvgProSolCB(i,j),stdb
               call printLine(RESULTS,printstring)
               write(printstring,21) AvgProSolCF(i,j)-AvgProSolC(i,j),AvgProSolCB(i,j)-AvgProSolC(i,j)
               call printLine(RESULTS,printstring)

c              lj
               std = stdDev(AvgProSolLJ2(i,j),AvgProSolLJ(i,j))
               stdf = stdDev(AvgProSolLJF2(i,j),AvgProSolLJF(i,j))
               stdb = stdDev(AvgProSolLJB2(i,j),AvgProSolLJB(i,j))
               write(printstring,20) 'LJ ',AvgProSolLJ(i,j),std,AvgProSolLJF(i,j),stdf,
     .                                                      AvgProSolLJB(i,j),stdb
               call printLine(RESULTS,printstring)
               write(printstring,21) AvgProSolLJF(i,j)-AvgProSolLJ(i,j),AvgProSolLJB(i,j)-AvgProSolLJ(i,j)
               call printLine(RESULTS,printstring)
             endif
           enddo
         enddo
       endif
c"       
       if (NProteins.gt.0 .and. NSolvents.gt.0) then
         do i=1,NProteins
           lgth = lineLength(ProteinName(i))
           write(printstring,33) i,ProteinName(i)(1:lgth)
33         format('Average protein-solvent energy for ',i4,' ',a)      
           call prettyLine(RESULTS,printstring)

c            coulomb
c             PRINT*, "en" ,AvgProSvnC
             std = stdDev(AvgProSvnC2(i),AvgProSvnC(i))
             stdf = stdDev(AvgProSvnCF2(i),AvgProSvnCF(i))
             stdb = stdDev(AvgProSvnCB2(i),AvgProSvnCB(i))
             write(printstring,20) 'COU',AvgProSvnC(i),std,AvgProSvnCF(i),stdf,
     .                                                      AvgProSvnCB(i),stdb
             call printLine(RESULTS,printstring)
             write(printstring,21) AvgProSvnCF(i)-AvgProSvnC(i),AvgProSvnCB(i)-AvgProSvnC(i)
             call printLine(RESULTS,printstring)

c            lj
             std = stdDev(AvgProSvnLJ2(i),AvgProSvnLJ(i))
             stdf = stdDev(AvgProSvnLJF2(i),AvgProSvnLJF(i))
             stdb = stdDev(AvgProSvnLJB2(i),AvgProSvnLJB(i))
             write(printstring,20) 'LJ ',AvgProSvnLJ(i),std,AvgProSvnLJF(i),stdf,
     .                                                    AvgProSvnLJB(i),stdb
             call printLine(RESULTS,printstring)
             write(printstring,21) AvgProSvnLJF(i)-AvgProSvnLJ(i),AvgProSvnLJB(i)-AvgProSvnLJ(i)
             call printLine(RESULTS,printstring)
         enddo
       endif

c        PRINT out the GC-solute stuff
c        if (NProteins.gt.0 .and. NGCSolutes.gt.0) then
c         do i=1,NProteins
c           lgth = lineLength(ProteinName(i))
c           write(printstring,33) i,ProteinName(i)(1:lgth)
c33         format('Average protein-GCsolute energy for ',i4,' ',a)      
c           call prettyLine(RESULTS,printstring)

c            coulomb
c             PRINT*, "en" ,AvgProSvnC
c             std = stdDev(AvgProSvnC2(i),AvgProSvnC(i))
c             stdf = stdDev(AvgProSvnCF2(i),AvgProSvnCF(i))
c             stdb = stdDev(AvgProSvnCB2(i),AvgProSvnCB(i))
c             write(printstring,20) 'COU',AvgProSvnC(i),std,AvgProSvnCF(i),stdf,
c     .                                                      AvgProSvnCB(i),stdb
c             call printLine(RESULTS,printstring)
c             write(printstring,21) AvgProSvnCF(i)-AvgProSvnC(i),AvgProSvnCB(i)-AvgProSvnC(i)
c             call printLine(RESULTS,printstring)
c
c            lj
c             std = stdDev(AvgProSvnLJ2(i),AvgProSvnLJ(i))
c             stdf = stdDev(AvgProSvnLJF2(i),AvgProSvnLJF(i))
c             stdb = stdDev(AvgProSvnLJB2(i),AvgProSvnLJB(i))
c             write(printstring,20) 'LJ',AvgProSvnLJ(i),std,AvgProSvnLJF(i),stdf,
c     .                                                    AvgProSvnLJB(i),stdb
c             call printLine(RESULTS,printstring)
c             write(printstring,21) AvgProSvnLJF(i)-AvgProSvnLJ(i),AvgProSvnLJB(i)-AvgProSvnLJ(i)
c             call printLine(RESULTS,printstring)
c         enddo
c       endif

       
       if (NSolutes.gt.1) then
         do i=1,NSolutes-1
           if (.not.IgnoreSolute(i)) then
             do j=i+1,NSolutes
               if (.not.IgnoreSolute(j)) then
                 lgth = lineLength(SoluteName(i))
                 lgth2 = lineLength(SoluteName(j))
                 write(printstring,34) i,SoluteName(i)(1:lgth),j,SoluteName(j)(1:lgth2)
34               format('Average solute-solute energy for ',i4,' ',a,' - ',i4,' ',a)             
                 call prettyLine(RESULTS,printstring)

c                coulomb
                 std = stdDev(AvgSolSolC2(i,j),AvgSolSolC(i,j))
                 stdf = stdDev(AvgSolSolCF2(i,j),AvgSolSolCF(i,j))
                 stdb = stdDev(AvgSolSolCB2(i,j),AvgSolSolCB(i,j))
                 write(printstring,20) 'COU',AvgSolSolC(i,j),std,AvgSolSolCF(i,j),stdf,
     .                                                           AvgSolSolCB(i,j),stdb
                 call printLine(RESULTS,printstring)
                 write(printstring,21) AvgSolSolCF(i,j)-AvgSolSolC(i,j),AvgSolSolCB(i,j)-AvgSolSolC(i,j)
                 call printLine(RESULTS,printstring)

c                lj
                 std = stdDev(AvgSolSolLJ2(i,j),AvgSolSolLJ(i,j))
                 stdf = stdDev(AvgSolSolLJF2(i,j),AvgSolSolLJF(i,j))
                 stdb = stdDev(AvgSolSolLJB2(i,j),AvgSolSolLJB(i,j))
                 write(printstring,20) 'LJ ',AvgSolSolLJ(i,j),std,AvgSolSolLJF(i,j),stdf,
     .                                                      AvgSolSolLJB(i,j),stdb
                 call printLine(RESULTS,printstring)
                 write(printstring,21) AvgSolSolLJF(i,j)-AvgSolSolLJ(i,j),AvgSolSolLJB(i,j)-AvgSolSolLJ(i,j)
                 call printLine(RESULTS,printstring)
               endif
             enddo
           endif
         enddo
       endif
      
       if (NSolutes.gt.0 .and. NSolvents.gt.0) then
         do i=1,NSolutes
           if (.not.IgnoreSolute(i)) then
             lgth = lineLength(SoluteName(i))
             write(printstring,35) i,SoluteName(i)(1:lgth)
35           format('Average solute-solvent energy for ',i4,' ',a20)
             call prettyLine(RESULTS,printstring)

c           coulomb
            std = stdDev(AvgSolSvnC2(i),AvgSolSvnC(i))
            stdf = stdDev(AvgSolSvnCF2(i),AvgSolSvnCF(i))
            stdb = stdDev(AvgSolSvnCB2(i),AvgSolSvnCB(i))
            write(printstring,20) 'COU',AvgSolSvnC(i),std,AvgSolSvnCF(i),stdf,
     .                                                    AvgSolSvnCB(i),stdb
            call printLine(RESULTS,printstring)
            write(printstring,21) AvgSolSvnCF(i)-AvgSolSvnC(i),AvgSolSvnCB(i)-AvgSolSvnC(i)
            call printLine(RESULTS,printstring)

c           lj
            std = stdDev(AvgSolSvnLJ2(i),AvgSolSvnLJ(i))
            stdf = stdDev(AvgSolSvnLJF2(i),AvgSolSvnLJF(i))
            stdb = stdDev(AvgSolSvnLJB2(i),AvgSolSvnLJB(i))
            write(printstring,20) 'LJ ',AvgSolSvnLJ(i),std,AvgSolSvnLJF(i),stdf,
     .                                                     AvgSolSvnLJB(i),stdb
            call printLine(RESULTS,printstring)
            write(printstring,21) AvgSolSvnLJF(i)-AvgSolSvnLJ(i),AvgSolSvnLJB(i)-AvgSolSvnLJ(i)
            call printLine(RESULTS,printstring)
           endif
         enddo
       endif
       
       if (NSolutes.gt.0 .and. NGCSolutes.gt.0) then
         do i=1,NSolutes
           if (.not.IgnoreSolute(i)) then
             lgth = lineLength(SoluteName(i))
             write(printstring,335) i,SoluteName(i)(1:lgth)
335          format('Average solute-GCS energy for ',i4,' ',a20)
             call prettyLine(RESULTS,printstring)

c           coulomb
            std = stdDev(AvgSolGCSC2(i),AvgSolGCSC(i))
            stdf = stdDev(AvgSolGCSCF2(i),AvgSolGCSCF(i))
            stdb = stdDev(AvgSolGCSCB2(i),AvgSolGCSCB(i))
            write(printstring,20) 'COU',AvgSolGCSC(i),std,AvgSolGCSCF(i),stdf,
     .                                                    AvgSolGCSCB(i),stdb
            call printLine(RESULTS,printstring)
            write(printstring,21) AvgSolGCSCF(i)-AvgSolGCSC(i),AvgSolGCSCB(i)-AvgSolGCSC(i)
            call printLine(RESULTS,printstring)

c           lj
            std = stdDev(AvgSolGCSLJ2(i),AvgSolGCSLJ(i))
            stdf = stdDev(AvgSolGCSLJF2(i),AvgSolGCSLJF(i))
            stdb = stdDev(AvgSolGCSLJB2(i),AvgSolGCSLJB(i))
            write(printstring,20) 'LJ ',AvgSolGCSLJ(i),std,AvgSolGCSLJF(i),stdf,
     .                                                     AvgSolGCSLJB(i),stdb
            call printLine(RESULTS,printstring)
            write(printstring,21) AvgSolGCSLJF(i)-AvgSolGCSLJ(i),AvgSolGCSLJB(i)-AvgSolGCSLJ(i)
            call printLine(RESULTS,printstring)
           endif
         enddo
       endif

c      print out the averages of the extra energy terms
      call printAverageExtraEnergies

c     print out free energy stuff
      call dGEnergies
 
      if (NThetaSolutes.gt.0) then
        write(printstring,"(A,I3,A)") " Individual theta solute values for ",NThetaSolutes," solutes"
        call prettyLine(RESULTS,printstring)
        do i=1,NThetaSolutes
          write(printstring,"('theta',i8,1x,f8.6)") i,SoluteTheta(i)
          call printLine(RESULTS,printstring)
        enddo
      endif
      
       if (NGCSolutes.gt.0) then
        write(printstring,*) "Individual theta values"
       call prettyLine(RESULTS,printstring)
        endif


      do i=1,NGCSolutes
       write(printstring,"('theta',i8,1x,f8.6)") i,thetaparam(i)
       call printLine(RESULTS,printstring)
      enddo

c      the variable 'space' determines 'cavCoordDepth'. By default both are zero
       if (NGCSolutes.gt.0.and.cavCoordDepth.gt.0) then
        write(printstring,"(A,F4.2)") "JAWS-I sites with theta greater than ",thres
       call prettyLine(RESULTS,printstring)
        endif
       do X = 1,cavCoordDepth
c        print*, "savedCoords?" ,stackGridX(X,1),stackGridY(X,2),stackGridZ(X,3)
        write(printstring,"('THETA',i8,3f12.3)") X,stackGridX(X,1),stackGridY(X,2),stackGridZ(X,3)
        call printLine(RESULTS,printstring)
       enddo
       cavCoordDepth = 0

       if (NGCSolutes.gt.0) then
        write(printstring,"(A,F4.2,A)") "JAWS-2 moves with theta above ",thres,
     .                                  " - centre position"
       call prettyLine(RESULTS,printstring)
        endif

       if (everyjaws2) then
         initial = avHighstack+1
       else
         initial = highstack
       endif
       do s = initial,highstack
c        print*, "savedCoords?" ,stackGridX(X,1),stackGridY(X,2),stackGridZ(X,3)
        write(printstring,"('HIGH',i8,3f12.3)") s,stackGridX(s,1),stackGridY(s,2),stackGridZ(s,3)
        call printLine(RESULTS,printstring)
       enddo
       avHighstack = highstack

       if (NGCSolutes.gt.0) then
        write(printstring,"(A,F4.2,A)") "JAWS-2 moves with theta below ",1-thres,
     .                                  " - centre position"
       call prettyLine(RESULTS,printstring)
        endif

       if (everyjaws2) then
         initial = avLowstack+1
       else
         initial = lowstack
       endif
       do t = initial,lowstack
c        prnt*, "savedCoords?" ,stackGridX(X,1),stackGridY(X,2),stackGridZ(X,3)
        write(printstring,"('LOW ',i8,3f12.3)") t,stackGridXL(t,1),stackGridYL(t,2),stackGridZL(t,3)
        call printLine(RESULTS,printstring)
       enddo
       avLowstack = lowstack
       return
       end
