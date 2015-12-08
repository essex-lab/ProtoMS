
      subroutine printAllEnergies
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
c     This routine prints the instantaneous energies of the system
c      in contrast to printAverages that print the running averages
c
c     However, retain the format of the results file, all stds
c     are printed as zeroes.
c
c     (C) Samuel Genheden      August 2014
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


       write(printstring,"('Number of data steps',i8)") 1
       call printLine(RESULTS,printstring)

       if (doLambdaRE) then
         write(printstring,"('Lambda replica = ',i8)") mylamid+1
         call printLine(RESULTS,printstring)
       endif

       if (doTemperatureRE) then
         write(printstring,"('Temperature replica = ',i8)") mytemid+1
         call printLine(RESULTS,printstring)
       endif
       
       if (doTemperatureRE.and.doLambdaRE) then
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
       
         write(printstring,"('Average volume = ',f12.4,' (',f12.4,')')") AvgVolume,ZERO
         call printLine(RESULTS,printstring)
       endif

c      write out the random number seed used to generate this trajectory
       write(printstring,"('Random number seed = ',i12)") ranseed
       call printLine(RESULTS,printstring)
                                  
       if (DifferentLambdas) then
         gfor = OLDTOTALENERGYF-OLDTOTALENERGY
         gback = OLDTOTALENERGYB-OLDTOTALENERGY

         if (abs(gfor).lt.TINY) gfor = ZERO
         if (abs(gback).lt.TINY) gback = ZERO
       
         write(printstring,19) "Backwards",gback,ZERO
         call printLine(RESULTS,printstring)
         
         write(printstring,19) "Forwards ",gfor,ZERO
         call printLine(RESULTS,printstring)

19       format(a9,' Free Energy = ',f15.4,' ( ',f8.3,' )')
       endif
       write(printstring,*) "Averages. Stddev in brackets. Average | forwards backwards"
       call prettyLine(RESULTS,printstring)
       
       write(printstring,43) OLDTOTALENERGY,ZERO,OLDTOTALENERGYF,ZERO,OLDTOTALENERGYB,ZERO
43     format('Total Energy ',f15.4,' ( ',f8.3,' ) | ',f15.4,' ( ',f8.3,' ) ',f15.4,' ( ',f8.3,' )')
       call printLine(RESULTS,printstring)

       if (NProteins.gt.0) then
         do i=1,NProteins
           if (NBonds(i).gt.0) then
c           internal bond energy
             lgth = lineLength(ProteinName(i))
             write(printstring,*) "Internal energy averages for protein ",i," ",ProteinName(i)(1:lgth)
             call prettyLine(RESULTS,printstring)
             write(printstring,20) "BND",OldProBnd(i),ZERO,OldProBndF(i),ZERO,OldProBndB(i),ZERO
20           format(a3,1x,f15.4,' ( ',f8.3,' ) | ',f15.4,' ( ',f8.3,' ) ',f15.4,' ( ',f8.3,' )')
             call printLine(RESULTS,printstring)
             write(printstring,21) OldProBndF(i)-OldProBnd(i),OldProBndB(i)-OldProBnd(i)
             call printLine(RESULTS,printstring)
21           format('DIFF',31x,f15.4,14x,f15.4)
           endif
           
           if (NAngles(i).gt.0) then
c            internal angle energy
             write(printstring,20) "ANG",OldProAng(i),ZERO,OldProAngF(i),ZERO,OldProAngB(i),ZERO
             call printLine(RESULTS,printstring)
             write(printstring,21) OldProAngF(i)-OldProAng(i),OldProAngB(i)-OldProAng(i)
             call printLine(RESULTS,printstring)
           endif
                
           if (NDihedrals(i).gt.0) then               
c            internal dihedral energy
             write(printstring,20) "DIH",OldProDih(i),ZERO,OldProDihF(i),ZERO,OldProDihB(i),ZERO
             call printLine(RESULTS,printstring)
             write(printstring,21) OldProDihF(i)-OldProDih(i),OldProDihB(i)-OldProDih(i)
             call printLine(RESULTS,printstring)
           endif

           if (NUreyBradleys(i).gt.0) then               
c            internal Urey-Bradley energy
             write(printstring,20) "UBY",OldProUBy(i),ZERO,OldProUByF(i),ZERO,OldProUByB(i),ZERO
             call printLine(RESULTS,printstring)
             write(printstring,21) OldProUByF(i)-OldProUBy(i),OldProUByB(i)-OldProUBy(i)
             call printLine(RESULTS,printstring)
           endif

c          internal nbc energy
           write(printstring,20) "NBC",OldProNBC(i),ZERO,OldProNBCF(i),ZERO,OldProNBCB(i),ZERO
           call printLine(RESULTS,printstring)
           write(printstring,21) OldProNBCF(i)-OldProNBC(i),OldProNBCB(i)-OldProNBC(i)
           call printLine(RESULTS,printstring)

c          internal nblj energy
           write(printstring,20) "NBL",OldProNBLJ(i),ZERO,OldProNBLJF(i),ZERO,OldProNBLJB(i),ZERO
           call printLine(RESULTS,printstring)
           write(printstring,21) OldProNBLJF(i)-OldProNBLJ(i),OldProNBLJB(i)-OldProNBLJ(i)
           call printLine(RESULTS,printstring)
         enddo
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
               write(printstring,20) 'BND',OldSolBnd(i),ZERO,OldSolBndF(i),ZERO,OldSolBndB(i),ZERO
               call printLine(RESULTS,printstring)
               write(printstring,21) OldSolBndF(i)-OldSolBnd(i),OldSolBndB(i)-OldSolBnd(i)
               call printLine(RESULTS,printstring)
             endif
           
             if (NSolAngles(i).gt.0) then
c              angle
               write(printstring,20) 'ANG',OldSolAng(i),ZERO,OldSolAngF(i),ZERO,OldSolAngB(i),ZERO
               call printLine(RESULTS,printstring)
               write(printstring,21) OldSolAngF(i)-OldSolAng(i),OldSolAngB(i)-OldSolAng(i)
               call printLine(RESULTS,printstring)
             endif

             if (NSolDihedrals(i).gt.0) then
c              dihedral
               write(printstring,20) 'DIH',OldSolDih(i),ZERO,OldSolDihF(i),ZERO,OldSolDihB(i),ZERO
               call printLine(RESULTS,printstring)
               write(printstring,21) OldSolDihF(i)-OldSolDih(i),OldSolDihB(i)-OldSolDih(i)
               call printLine(RESULTS,printstring)
             endif

             if (NSolUreyBradleys(i).gt.0) then
c              Urey-Bradley
               write(printstring,20) 'UBY',OldSolUBy(i),ZERO,OldSolUByF(i),ZERO,OldSolUByB(i),ZERO
               call printLine(RESULTS,printstring)
               write(printstring,21) OldSolUByF(i)-OldSolUBy(i),OldSolUByB(i)-OldSolUBy(i)
               call printLine(RESULTS,printstring)
             endif

c            nbc
             write(printstring,20) 'NBC',OldSolNBC(i),ZERO,OldSolNBCF(i),ZERO,OldSolNBCB(i),ZERO
             call printLine(RESULTS,printstring)
             write(printstring,21) OldSolNBCF(i)-OldSolNBC(i),OldSolNBCB(i)-OldSolNBC(i)
             call printLine(RESULTS,printstring)

c            nblj
             write(printstring,20) 'NBL',OldSolNBLJ(i),ZERO,OldSolNBLJF(i),ZERO,OldSolNBLJB(i),ZERO
             call printLine(RESULTS,printstring)
             write(printstring,21) OldSolNBLJF(i)-OldSolNBLJ(i),OldSolNBLJB(i)-OldSolNBLJ(i)
             call printLine(RESULTS,printstring)
           endif
         enddo
       endif
       
       if (NSolvents.gt.0) then
         call prettyLine(RESULTS,"Average solvent-solvent energy")
c        coulomb energy         
         write(printstring,20) 'COU',OldSvnSvnC,ZERO,OldSvnSvnCF,ZERO,OldSvnSvnCB,ZERO
         call printLine(RESULTS,printstring)
         write(printstring,21) OldSvnSvnCF-OldSvnSvnC,OldSvnSvnCB-OldSvnSvnC
         call printLine(RESULTS,printstring)
         
c        lj energy
         write(printstring,20) 'LJ ',OldSvnSvnLJ,ZERO,OldSvnSvnLJF,ZERO,OldSvnSvnLJB,ZERO
         call printLine(RESULTS,printstring)
         write(printstring,21) OldSvnSvnLJF-OldSvnSvnLJ,OldSvnSvnLJB-OldSvnSvnLJ
         call printLine(RESULTS,printstring)
       endif
       
       if (NGCSolutes.gt.0) then
         call prettyLine(RESULTS,"Average GCS-GCS energy")
c        coulomb energy         
         write(printstring,20) 'COU',OldGCSGCSC,ZERO,OldGCSGCSCF,ZERO,OldGCSGCSCB,ZERO
         call printLine(RESULTS,printstring)
         write(printstring,21) OldGCSGCSCF-OldGCSGCSC,OldGCSGCSCB-OldGCSGCSC
         call printLine(RESULTS,printstring)
         
c        lj energy
         write(printstring,20) 'LJ ',OldGCSGCSLJ,ZERO,oldGCSGCSLJF,Zero,OldGCSGCSLJB,ZERO
         call printLine(RESULTS,printstring)
         write(printstring,21) OldGCSGCSLJF-OldGCSGCSLJ,OldGCSGCSLJB-OldGCSGCSLJ
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
             write(printstring,20) 'COU',OLDProGCSC(i),ZERO,OldProGCSCF(i),ZERO,
     .                                                      OldProGCSCB(i),ZERO
             call printLine(RESULTS,printstring)
             write(printstring,21) OldProGCSCF(i)-OldProGCSC(i),OldProGCSCB(i)-OldProGCSC(i)
             call printLine(RESULTS,printstring)

c            lj
             write(printstring,20) 'LJ',OldProGCSLJ(i),ZERO,OldProGCSLJF(i),ZERO,
     .                                                    OldProGCSLJB(i),ZERO
             call printLine(RESULTS,printstring)
             write(printstring,21) OldProGCSLJF(i)-OldProGCSLJ(i),OldProGCSLJB(i)-OldProGCSLJ(i)
             call printLine(RESULTS,printstring)
         enddo
       endif

cc      if (NGCSolutes.gt.0 .and. NSolvents.gt.0) then
c         do i=1,NGCSolutes
c           lgth = lineLength(GCSName(i))
c           write(printstring,932) i,GCSName(i)(1:lgth)
c932        format('Average GCS-solvent energy for ',i4,' ',a)      
c           call prettyLine(RESULTS,printstring)
c
cc            coulomb
cc             PRINT*, "en" ,AvgGCSSvnC
c             std = stdDev(AvgGCSSvnC2(i),AvgGCSSvnC(i))
c             stdf = stdDev(AvgGCSSvnCF2(i),AvgGCSSvnCF(i))
c             stdb = stdDev(AvgGCSSvnCB2(i),AvgGCSSvnCB(i))
c             write(printstring,20) 'COU',AvgGCSSvnC(i),std,AvgGCSSvnCF(i),stdf,
c     .                                                      AvgGCSSvnCB(i),stdb
c             call printLine(RESULTS,printstring)
c             write(printstring,21) AvgGCSSvnCF(i)-AvgGCSSvnC(i),AvgGCSSvnCB(i)-AvgGCSSvnC(i)
c             call printLine(RESULTS,printstring)
c
cc            lj
c             std = stdDev(AvgGCSSvnLJ2(i),AvgGCSSvnLJ(i))
c             stdf = stdDev(AvgGCSSvnLJF2(i),AvgGCSSvnLJF(i))
c             stdb = stdDev(AvgGCSSvnLJB2(i),AvgGCSSvnLJB(i))
c             write(printstring,20) 'LJ',AvgGCSSvnLJ(i),std,AvgGCSSvnLJF(i),stdf,
c     .                                                    AvgGCSSvnLJB(i),stdb
c             call printLine(RESULTS,printstring)
c             write(printstring,21) AvgGCSSvnLJF(i)-AvgGCSSvnLJ(i),AvgGCSSvnLJB(i)-AvgGCSSvnLJ(i)
c             call printLine(RESULTS,printstring)
c         enddo
c       endif

       
          
c"

        

         do i=1,NProteins-1
           do j=i+1,NProteins
             write(printstring,30) i,ProteinName(i)(1:20),j,ProteinName(j)(1:20)
30           format('Average protein-protein energy for ',i4,' ',a20,' - ',i4,' ',a20)             
             call prettyLine(RESULTS,printstring)

c            coulomb
             write(printstring,20) 'COU',OldProProC(i,j),ZERO,OldProProCF(i,j),ZERO,
     .                                                      OldProProCB(i,j),ZERO
             call printLine(RESULTS,printstring)
             write(printstring,21) OldProProCF(i,j)-OldProProC(i,j),OldProProCB(i,j)-OldProProC(i,j)
             call printLine(RESULTS,printstring)

c            lj
             write(printstring,20) 'LJ ',OldProProLJ(i,j),ZERO,OldProProLJF(i,j),ZERO,
     .                                                      OldProProLJB(i,j),ZERO
             call printLine(RESULTS,printstring)
             write(printstring,21) OldProProLJF(i,j)-OldProProLJ(i,j),OldProProLJB(i,j)-OldProProLJ(i,j)
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
               write(printstring,20) 'COU',OldProSolC(i,j),ZERO,OldProSolCF(i,j),ZERO,
     .                                                      OldProSolCB(i,j),ZERO
               call printLine(RESULTS,printstring)
               write(printstring,21) OldProSolCF(i,j)-OldProSolC(i,j),OldProSolCB(i,j)-OldProSolC(i,j)
               call printLine(RESULTS,printstring)

c              lj
               write(printstring,20) 'LJ ',OldProSolLJ(i,j),ZERO,OldProSolLJF(i,j),ZERO,
     .                                                      OldProSolLJB(i,j),ZERO
               call printLine(RESULTS,printstring)
               write(printstring,21) OldProSolLJF(i,j)-OldProSolLJ(i,j),OldProSolLJB(i,j)-OldProSolLJ(i,j)
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
             write(printstring,20) 'COU',OldProSvnC(i),ZERO,OldProSvnCF(i),ZERO,
     .                                                      OldProSvnCB(i),ZERO
             call printLine(RESULTS,printstring)
             write(printstring,21) OldProSvnCF(i)-OldProSvnC(i),OldProSvnCB(i)-OldProSvnC(i)
             call printLine(RESULTS,printstring)

c            lj
             write(printstring,20) 'LJ',OldProSvnLJ(i),ZERO,OldProSvnLJF(i),ZERO,
     .                                                    OldProSvnLJB(i),ZERO
             call printLine(RESULTS,printstring)
             write(printstring,21) OldProSvnLJF(i)-OldProSvnLJ(i),OldProSvnLJB(i)-OldProSvnLJ(i)
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
                 write(printstring,20) 'COU',OldSolSolC(i,j),ZERO,OldSolSolCF(i,j),ZERO,
     .                                                           OldSolSolCB(i,j),ZERO
                 call printLine(RESULTS,printstring)
                 write(printstring,21) OldSolSolCF(i,j)-OldSolSolC(i,j),OldSolSolCB(i,j)-OldSolSolC(i,j)
                 call printLine(RESULTS,printstring)

c                lj
                 write(printstring,20) 'LJ ',OldSolSolLJ(i,j),ZERO,OldSolSolLJF(i,j),ZERO,
     .                                                      OldSolSolLJB(i,j),ZERO
                 call printLine(RESULTS,printstring)
                 write(printstring,21) OldSolSolLJF(i,j)-OldSolSolLJ(i,j),OldSolSolLJB(i,j)-OldSolSolLJ(i,j)
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
            write(printstring,20) 'COU',OldSolSvnC(i),ZERO,OldSolSvnCF(i),ZERO,
     .                                                    OldSolSvnCB(i),ZERO
            call printLine(RESULTS,printstring)
            write(printstring,21) OldSolSvnCF(i)-OldSolSvnC(i),OldSolSvnCB(i)-OldSolSvnC(i)
            call printLine(RESULTS,printstring)

c           lj
            write(printstring,20) 'LJ ',OldSolSvnLJ(i),ZERO,OldSolSvnLJF(i),ZERO,
     .                                                     OldSolSvnLJB(i),ZERO
            call printLine(RESULTS,printstring)
            write(printstring,21) OldSolSvnLJF(i)-OldSolSvnLJ(i),OldSolSvnLJB(i)-OldSolSvnLJ(i)
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
            write(printstring,20) 'COU',OldSolGCSC(i),ZERO,OldSolGCSCF(i),ZERO,
     .                                                    OldSolGCSCB(i),ZERO
            call printLine(RESULTS,printstring)
            write(printstring,21) OldSolGCSCF(i)-OldSolGCSC(i),OldSolGCSCB(i)-OldSolGCSC(i)
            call printLine(RESULTS,printstring)

c           lj
            write(printstring,20) 'LJ ',OldSolGCSLJ(i),ZERO,OldSolGCSLJF(i),ZERO,
     .                                                     OldSolGCSLJB(i),ZERO
            call printLine(RESULTS,printstring)
            write(printstring,21) OldSolGCSLJF(i)-OldSolGCSLJ(i),OldSolGCSLJB(i)-OldSolGCSLJ(i)
            call printLine(RESULTS,printstring)
           endif
         enddo
       endif
       
c      print out the averages of the extra energy terms
       call printAllExtraEnergies
 
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
         initial = instHighstack+1
       else
         initial = highstack
       endif
       do s = initial,highstack
c        print*, "savedCoords?" ,stackGridX(X,1),stackGridY(X,2),stackGridZ(X,3)
        write(printstring,"('HIGH',i8,3f12.3)") s,stackGridX(s,1),stackGridY(s,2),stackGridZ(s,3)
        call printLine(RESULTS,printstring)
       enddo
       instHighstack = highstack

       if (NGCSolutes.gt.0) then
        write(printstring,"(A,F4.2,A)") "JAWS-2 moves with theta below ",1-thres,
     .                                  " - centre position"
       call prettyLine(RESULTS,printstring)
        endif

       if (everyjaws2) then
         initial = instLowstack+1
       else
         initial = lowstack
       endif
       do t = initial,lowstack
c        prnt*, "savedCoords?" ,stackGridX(X,1),stackGridY(X,2),stackGridZ(X,3)
        write(printstring,"('LOW ',i8,3f12.3)") t,stackGridXL(t,1),stackGridYL(t,2),stackGridZL(t,3)
        call printLine(RESULTS,printstring)
       enddo
       instLowstack = lowstack
       return
       end
