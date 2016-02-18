      
      subroutine calcAverageEnergy
      include 'dimensions.inc'
      include 'metropolis.inc'
      include 'averages.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'solvents.inc'
      include 'flags.inc'
      include 'energies.inc'
      include 'pressure.inc'
      include 'extraenergies.inc'
      include 'printer.inc'
c#####################################################
c
c     This calculates averages of important
c     quantities - especially free energies!
c       Updated MSB 27/04/09
c#####################################################
      
      double precision bigratio,smallratio,rsteps
      double precision gtemp      
      integer i,j
                                            
c     averages calculated as a running value, based on the equation
c     
c           newAverage = (oldAverage * nsteps-1/nsteps) + newValue/nsteps
c
c       in this way, we avoid large running totals which have to be divided
c       (and could be subject to error), and we are always dealing with numbers
c       that should be well handled as doubles (since nsteps-1/nsteps will 
c       lie between 0.5 and 1, and 1/nsteps will lie between 0 and 0.5. 
c       Also, 1/nsteps = 1 - (nsteps-1/nsteps)

c     first, increment the number of steps in the average
      nsteps = nsteps + 1

c     calculate (nsteps-1/nsteps) and (1/nsteps)
      rsteps = dble(nsteps)
      bigratio = (rsteps-ONE) / rsteps
      smallratio = ONE - bigratio
c      write(printstring,*) 'Adding to average ',OldTotalEnergy
c      call printLine(DEBUG,printstring)
c     calculate the average total energy (from OLDTOTALENERGY)

c      write(printstring,*) 'Total Energy is ',OldTotalEnergy
c      call printLine(DEBUG,printstring)

      AVGTOTALENERGY = (AVGTOTALENERGY*bigratio) + (OLDTOTALENERGY*smallratio)
      AVGTOTALENERGYF = (AVGTOTALENERGYF*bigratio) + (OLDTOTALENERGYF*smallratio)
      AVGTOTALENERGYB = (AVGTOTALENERGYB*bigratio) + (OLDTOTALENERGYB*smallratio)
c      PRINT*, "OLDTOTALENERGY" ,OLDTOTALENERGY
c      PRINT*, "AVGTOTALENERGY" ,AVGTOTALENERGY
c      PRINT*, "bigratio" ,bigratio
c     also need the stddev, defined as 
c         stddev**2 = mean(squares) - square(mean)      
c       thus we also need to average the square of the values...
      AVGTOTALENERGY2 = (AVGTOTALENERGY2*bigratio) + (smallratio*(OLDTOTALENERGY**2))
      AVGTOTALENERGYF2 = (AVGTOTALENERGYF2*bigratio) + (smallratio*(OLDTOTALENERGYF**2))
      AVGTOTALENERGYB2 = (AVGTOTALENERGYB2*bigratio) + (smallratio*(OLDTOTALENERGYB**2))

c     average the volume
      AvgVolume = (AvgVolume*bigratio) + (smallratio*Volume)
      AvgVolume2 = (AvgVolume2*bigratio) + (smallratio*(Volume**2))
         
c     if we are doing FEP, then calculate the average forwards and backwards
c     energies (see the zwanzig equation...
      if (DifferentLambdas) then
       
c       dG(A->B) = -kT ln < exp(-(Eb -Ea) / kT) >a
        gtemp = dexp( -(OLDTOTALENERGYB-OLDTOTALENERGY) * mt_beta_T )
        AVGBACKENERGY = (AVGBACKENERGY*bigratio) + (gtemp*smallratio)
        
        gtemp = OLDTOTALENERGYB-OLDTOTALENERGY
        AVGBACKENERGY1 = (AVGBACKENERGY1*bigratio) + (smallratio*gtemp)
        AVGBACKENERGY2 = (AVGBACKENERGY2*bigratio) + (smallratio*(gtemp**2))
        
        gtemp = dexp( -(OLDTOTALENERGYF-OLDTOTALENERGY) * mt_beta_T )
        AVGFORWARDSENERGY = (AVGFORWARDSENERGY*bigratio) + (gtemp*smallratio)
        
        gtemp = OLDTOTALENERGYF-OLDTOTALENERGY
        AVGFORWARDSENERGY1 = (AVGFORWARDSENERGY1*bigratio) + (smallratio*gtemp)
        AVGFORWARDSENERGY2 = (AVGFORWARDSENERGY2*bigratio) + (smallratio*(gtemp**2))

      endif

c     now calculate the averages/stddev of the simulation energy components
c       first, all of the protein related averages      
      if (NProteins.gt.0) then
        do i=1,NProteins
        
          AvgProBnd(i) = (AvgProBnd(i)*bigratio) + (OldProBnd(i)*smallratio)
          AvgProAng(i) = (AvgProAng(i)*bigratio) + (OldProAng(i)*smallratio)
          AvgProDih(i) = (AvgProDih(i)*bigratio) + (OldProDih(i)*smallratio)
          AvgProUBy(i) = (AvgProUBy(i)*bigratio) + (OldProUBy(i)*smallratio)
          AvgProNBC(i) = (AvgProNBC(i)*bigratio) + (OldProNBC(i)*smallratio)
          AvgProNBLJ(i) = (AvgProNBLJ(i)*bigratio) + (OldProNBLJ(i)*smallratio)
          AvgProBnd2(i) = (AvgProBnd2(i)*bigratio) + (OldProBnd(i)**2 * smallratio)
          AvgProAng2(i) = (AvgProAng2(i)*bigratio) + (OldProAng(i)**2 * smallratio)
          AvgProDih2(i) = (AvgProDih2(i)*bigratio) + (OldProDih(i)**2 * smallratio)
          AvgProUBy2(i) = (AvgProUBy2(i)*bigratio) + (OldProUBy(i)**2 * smallratio)
          AvgProNBC2(i) = (AvgProNBC2(i)*bigratio) + (OldProNBC(i)**2 * smallratio)
          AvgProNBLJ2(i) = (AvgProNBLJ2(i)*bigratio) + (OldProNBLJ(i)**2 * smallratio)
          AvgProBndF(i) = (AvgProBndF(i)*bigratio) + (OldProBndF(i)*smallratio)
          AvgProAngF(i) = (AvgProAngF(i)*bigratio) + (OldProAngF(i)*smallratio)
          AvgProDihF(i) = (AvgProDihF(i)*bigratio) + (OldProDihF(i)*smallratio)
          AvgProUByF(i) = (AvgProUByF(i)*bigratio) + (OldProUByF(i)*smallratio)
          AvgProNBCF(i) = (AvgProNBCF(i)*bigratio) + (OldProNBCF(i)*smallratio)
          AvgProNBLJF(i) = (AvgProNBLJF(i)*bigratio) + (OldProNBLJF(i)*smallratio)
          AvgProBndF2(i) = (AvgProBndF2(i)*bigratio) + (OldProBndF(i)**2 * smallratio)
          AvgProAngF2(i) = (AvgProAngF2(i)*bigratio) + (OldProAngF(i)**2 * smallratio)
          AvgProDihF2(i) = (AvgProDihF2(i)*bigratio) + (OldProDihF(i)**2 * smallratio)
          AvgProUByF2(i) = (AvgProUByF2(i)*bigratio) + (OldProUByF(i)**2 * smallratio)
          AvgProNBCF2(i) = (AvgProNBCF2(i)*bigratio) + (OldProNBCF(i)**2 * smallratio)
          AvgProNBLJF2(i) = (AvgProNBLJF2(i)*bigratio) + (OldProNBLJF(i)**2 * smallratio)
          AvgProBndB(i) = (AvgProBndB(i)*bigratio) + (OldProBndB(i)*smallratio)
          AvgProAngB(i) = (AvgProAngB(i)*bigratio) + (OldProAngB(i)*smallratio)
          AvgProDihB(i) = (AvgProDihB(i)*bigratio) + (OldProDihB(i)*smallratio)
          AvgProUByB(i) = (AvgProUByB(i)*bigratio) + (OldProUByB(i)*smallratio)
          AvgProNBCB(i) = (AvgProNBCB(i)*bigratio) + (OldProNBCB(i)*smallratio)
          AvgProNBLJB(i) = (AvgProNBLJB(i)*bigratio) + (OldProNBLJB(i)*smallratio)
          AvgProBndB2(i) = (AvgProBndB2(i)*bigratio) + (OldProBndB(i)**2 * smallratio)
          AvgProAngB2(i) = (AvgProAngB2(i)*bigratio) + (OldProAngB(i)**2 * smallratio)
          AvgProDihB2(i) = (AvgProDihB2(i)*bigratio) + (OldProDihB(i)**2 * smallratio)
          AvgProUByB2(i) = (AvgProUByB2(i)*bigratio) + (OldProUByB(i)**2 * smallratio)
          AvgProNBCB2(i) = (AvgProNBCB2(i)*bigratio) + (OldProNBCB(i)**2 * smallratio)
          AvgProNBLJB2(i) = (AvgProNBLJB2(i)*bigratio) + (OldProNBLJB(i)**2 * smallratio)
          
          if (NProteins.gt.1) then
            do j=1,NProteins
              AvgProProC(i,j) = (AvgProProC(i,j)*bigratio) + (OldProProC(i,j)*smallratio)
              AvgProProLJ(i,j) = (AvgProProLJ(i,j)*bigratio) + (OldProProLJ(i,j)*smallratio)
              AvgProProC2(i,j) = (AvgProProC2(i,j)*bigratio) + (OldProProC(i,j)**2 * smallratio)
              AvgProProLJ2(i,j) = (AvgProProLJ2(i,j)*bigratio) + (OldProProLJ(i,j)**2 * smallratio)
              AvgProProCF(i,j) = (AvgProProCF(i,j)*bigratio) + (OldProProCF(i,j)*smallratio)
              AvgProProLJF(i,j) = (AvgProProLJF(i,j)*bigratio) + (OldProProLJF(i,j)*smallratio)
              AvgProProCF2(i,j) = (AvgProProCF2(i,j)*bigratio) + (OldProProCF(i,j)**2 * smallratio)
              AvgProProLJF2(i,j) = (AvgProProLJF2(i,j)*bigratio) + (OldProProLJF(i,j)**2 * smallratio)
              AvgProProCB(i,j) = (AvgProProCB(i,j)*bigratio) + (OldProProCB(i,j)*smallratio)
              AvgProProLJB(i,j) = (AvgProProLJB(i,j)*bigratio) + (OldProProLJB(i,j)*smallratio)
              AvgProProCB2(i,j) = (AvgProProCB2(i,j)*bigratio) + (OldProProCB(i,j)**2 * smallratio)
              AvgProProLJB2(i,j) = (AvgProProLJB2(i,j)*bigratio) + (OldProProLJB(i,j)**2 * smallratio)
            enddo
          endif
              
          if (NSolutes.gt.0) then
            do j=1,NSolutes
              AvgProSolC(i,j) = (AvgProSolC(i,j)*bigratio) + (OldProSolC(i,j)*smallratio)
              AvgProSolLJ(i,j) = (AvgProSolLJ(i,j)*bigratio) + (OldProSolLJ(i,j)*smallratio)
              AvgProSolC2(i,j) = (AvgProSolC2(i,j)*bigratio) + (OldProSolC(i,j)**2 * smallratio)
              AvgProSolLJ2(i,j) = (AvgProSolLJ2(i,j)*bigratio) + (OldProSolLJ(i,j)**2 * smallratio)
              AvgProSolCF(i,j) = (AvgProSolCF(i,j)*bigratio) + (OldProSolCF(i,j)*smallratio)
              AvgProSolLJF(i,j) = (AvgProSolLJF(i,j)*bigratio) + (OldProSolLJF(i,j)*smallratio)
              AvgProSolCF2(i,j) = (AvgProSolCF2(i,j)*bigratio) + (OldProSolCF(i,j)**2 * smallratio)
              AvgProSolLJF2(i,j) = (AvgProSolLJF2(i,j)*bigratio) + (OldProSolLJF(i,j)**2 * smallratio)
              AvgProSolCB(i,j) = (AvgProSolCB(i,j)*bigratio) + (OldProSolCB(i,j)*smallratio)
              AvgProSolLJB(i,j) = (AvgProSolLJB(i,j)*bigratio) + (OldProSolLJB(i,j)*smallratio)
              AvgProSolCB2(i,j) = (AvgProSolCB2(i,j)*bigratio) + (OldProSolCB(i,j)**2 * smallratio)
              AvgProSolLJB2(i,j) = (AvgProSolLJB2(i,j)*bigratio) + (OldProSolLJB(i,j)**2 * smallratio)
            enddo
          endif
          
          if (NSolvents.gt.0) then
            AvgProSvnC(i) = (AvgProSvnC(i)*bigratio) + (OldProSvnC(i)*smallratio)
            AvgProSvnLJ(i) = (AvgProSvnLJ(i)*bigratio) + (OldProSvnLJ(i)*smallratio)
            AvgProSvnC2(i) = (AvgProSvnC2(i)*bigratio) + (OldProSvnC(i)**2 * smallratio)
            AvgProSvnLJ2(i) = (AvgProSvnLJ2(i)*bigratio) + (OldProSvnLJ(i)**2 * smallratio)
            AvgProSvnCF(i) = (AvgProSvnCF(i)*bigratio) + (OldProSvnCF(i)*smallratio)
            AvgProSvnLJF(i) = (AvgProSvnLJF(i)*bigratio) + (OldProSvnLJF(i)*smallratio)
            AvgProSvnCF2(i) = (AvgProSvnCF2(i)*bigratio) + (OldProSvnCF(i)**2 * smallratio)
            AvgProSvnLJF2(i) = (AvgProSvnLJF2(i)*bigratio) + (OldProSvnLJF(i)**2 * smallratio)
            AvgProSvnCB(i) = (AvgProSvnCB(i)*bigratio) + (OldProSvnCB(i)*smallratio)
            AvgProSvnLJB(i) = (AvgProSvnLJB(i)*bigratio) + (OldProSvnLJB(i)*smallratio)
            AvgProSvnCB2(i) = (AvgProSvnCB2(i)*bigratio) + (OldProSvnCB(i)**2 * smallratio)
            AvgProSvnLJB2(i) = (AvgProSvnLJB2(i)*bigratio) + (OldProSvnLJB(i)**2 * smallratio)
          endif
        
         if (NGCSolutes.gt.0) then
            AvgProGCSC(i) = (AvgProGCSC(i)*bigratio) + (OldProGCSC(i)*smallratio)
            AvgProGCSLJ(i) = (AvgProGCSLJ(i)*bigratio) + (OldProGCSLJ(i)*smallratio)
            AvgProGCSC2(i) = (AvgProGCSC2(i)*bigratio) + (OldProGCSC(i)**2 * smallratio)
            AvgProGCSLJ2(i) = (AvgProGCSLJ2(i)*bigratio) + (OldProGCSLJ(i)**2 * smallratio)
            AvgProGCSCF(i) = (AvgProGCSCF(i)*bigratio) + (OldProGCSCF(i)*smallratio)
            AvgProGCSCF2(i) = (AvgProGCSCF2(i)*bigratio) + (OldProGCSCF(i)**2 * smallratio)
            AvgProGCSLJF(i) = (AvgProGCSLJF(i)*bigratio) + (OldProGCSLJF(i)*smallratio)
            AvgProGCSLJB(i) = (AvgProGCSLJB(i)*bigratio) + (OldProGCSLJB(i)*smallratio)
            AvgProGCSLJF2(i) = (AvgProGCSLJF2(i)*bigratio) + (OldProGCSLJF(i)**2 * smallratio)
            AvgProGCSCB(i) = (AvgProGCSCB(i)*bigratio) + (OldProGCSCB(i)*smallratio)
            AvgProGCSCB2(i) = (AvgProGCSCB2(i)*bigratio) + (OldProGCSCB(i)**2 * smallratio)
            AvgProGCSLJB2(i) = (AvgProGCSLJB2(i)*bigratio) + (OldProGCSLJB(i)**2 * smallratio)
          endif

        enddo
      endif        

c       now the solute related averages
      if (NSolutes.gt.0) then
        do i=1,NSolutes
          AvgSolBnd(i) = (AvgSolBnd(i)*bigratio) + (OldSolBnd(i)*smallratio)
          AvgSolAng(i) = (AvgSolAng(i)*bigratio) + (OldSolAng(i)*smallratio)
          AvgSolDih(i) = (AvgSolDih(i)*bigratio) + (OldSolDih(i)*smallratio)
          AvgSolUBy(i) = (AvgSolUBy(i)*bigratio) + (OldSolUBy(i)*smallratio)
          AvgSolNBC(i) = (AvgSolNBC(i)*bigratio) + (OldSolNBC(i)*smallratio)
          AvgSolNBLJ(i) = (AvgSolNBLJ(i)*bigratio) + (OldSolNBLJ(i)*smallratio)
          AvgSolBnd2(i) = (AvgSolBnd2(i)*bigratio) + (OldSolBnd(i)**2 * smallratio)
          AvgSolAng2(i) = (AvgSolAng2(i)*bigratio) + (OldSolAng(i)**2 * smallratio)
          AvgSolDih2(i) = (AvgSolDih2(i)*bigratio) + (OldSolDih(i)**2 * smallratio)
          AvgSolUBy2(i) = (AvgSolUBy2(i)*bigratio) + (OldSolUBy(i)**2 * smallratio)
          AvgSolNBC2(i) = (AvgSolNBC2(i)*bigratio) + (OldSolNBC(i)**2 * smallratio)
          AvgSolNBLJ2(i) = (AvgSolNBLJ2(i)*bigratio) + (OldSolNBLJ(i)**2 * smallratio)
          AvgSolBndF(i) = (AvgSolBndF(i)*bigratio) + (OldSolBndF(i)*smallratio)
          AvgSolAngF(i) = (AvgSolAngF(i)*bigratio) + (OldSolAngF(i)*smallratio)
          AvgSolDihF(i) = (AvgSolDihF(i)*bigratio) + (OldSolDihF(i)*smallratio)
          AvgSolUByF(i) = (AvgSolUByF(i)*bigratio) + (OldSolUByF(i)*smallratio)
          AvgSolNBCF(i) = (AvgSolNBCF(i)*bigratio) + (OldSolNBCF(i)*smallratio)
          AvgSolNBLJF(i) = (AvgSolNBLJF(i)*bigratio) + (OldSolNBLJF(i)*smallratio)
          AvgSolBndF2(i) = (AvgSolBndF2(i)*bigratio) + (OldSolBndF(i)**2 * smallratio)
          AvgSolAngF2(i) = (AvgSolAngF2(i)*bigratio) + (OldSolAngF(i)**2 * smallratio)
          AvgSolDihF2(i) = (AvgSolDihF2(i)*bigratio) + (OldSolDihF(i)**2 * smallratio)
          AvgSolUByF2(i) = (AvgSolUByF2(i)*bigratio) + (OldSolUByF(i)**2 * smallratio)
          AvgSolNBCF2(i) = (AvgSolNBCF2(i)*bigratio) + (OldSolNBCF(i)**2 * smallratio)
          AvgSolNBLJF2(i) = (AvgSolNBLJF2(i)*bigratio) + (OldSolNBLJF(i)**2 * smallratio)
          AvgSolBndB(i) = (AvgSolBndB(i)*bigratio) + (OldSolBndB(i)*smallratio)
          AvgSolAngB(i) = (AvgSolAngB(i)*bigratio) + (OldSolAngB(i)*smallratio)
          AvgSolDihB(i) = (AvgSolDihB(i)*bigratio) + (OldSolDihB(i)*smallratio)
          AvgSolUByB(i) = (AvgSolUByB(i)*bigratio) + (OldSolUByB(i)*smallratio)
          AvgSolNBCB(i) = (AvgSolNBCB(i)*bigratio) + (OldSolNBCB(i)*smallratio)
          AvgSolNBLJB(i) = (AvgSolNBLJB(i)*bigratio) + (OldSolNBLJB(i)*smallratio)
          AvgSolBndB2(i) = (AvgSolBndB2(i)*bigratio) + (OldSolBndB(i)**2 * smallratio)
          AvgSolAngB2(i) = (AvgSolAngB2(i)*bigratio) + (OldSolAngB(i)**2 * smallratio)
          AvgSolDihB2(i) = (AvgSolDihB2(i)*bigratio) + (OldSolDihB(i)**2 * smallratio)
          AvgSolUByB2(i) = (AvgSolUByB2(i)*bigratio) + (OldSolUByB(i)**2 * smallratio)
          AvgSolNBCB2(i) = (AvgSolNBCB2(i)*bigratio) + (OldSolNBCB(i)**2 * smallratio)
          AvgSolNBLJB2(i) = (AvgSolNBLJB2(i)*bigratio) + (OldSolNBLJB(i)**2 * smallratio)

          if (NSolutes.gt.1) then
            do j=1,NSolutes
              AvgSolSolC(i,j) = (AvgSolSolC(i,j)*bigratio) + (OldSolSolC(i,j)*smallratio)
              AvgSolSolLJ(i,j) = (AvgSolSolLJ(i,j)*bigratio) + (OldSolSolLJ(i,j)*smallratio)
              AvgSolSolC2(i,j) = (AvgSolSolC2(i,j)*bigratio) + (OldSolSolC(i,j)**2*smallratio)
              AvgSolSolLJ2(i,j) = (AvgSolSolLJ2(i,j)*bigratio) + (OldSolSolLJ(i,j)**2*smallratio)
              AvgSolSolCF(i,j) = (AvgSolSolCF(i,j)*bigratio) + (OldSolSolCF(i,j)*smallratio)
              AvgSolSolLJF(i,j) = (AvgSolSolLJF(i,j)*bigratio) + (OldSolSolLJF(i,j)*smallratio)
              AvgSolSolCF2(i,j) = (AvgSolSolCF2(i,j)*bigratio) + (OldSolSolCF(i,j)**2*smallratio)
              AvgSolSolLJF2(i,j) = (AvgSolSolLJF2(i,j)*bigratio) + (OldSolSolLJF(i,j)**2*smallratio)
              AvgSolSolCB(i,j) = (AvgSolSolCB(i,j)*bigratio) + (OldSolSolCB(i,j)*smallratio)
              AvgSolSolLJB(i,j) = (AvgSolSolLJB(i,j)*bigratio) + (OldSolSolLJB(i,j)*smallratio)
              AvgSolSolCB2(i,j) = (AvgSolSolCB2(i,j)*bigratio) + (OldSolSolCB(i,j)**2*smallratio)
              AvgSolSolLJB2(i,j) = (AvgSolSolLJB2(i,j)*bigratio) + (OldSolSolLJB(i,j)**2*smallratio)
            enddo
          endif

          if (NSolvents.gt.0) then
            AvgSolSvnC(i) = (AvgSolSvnC(i)*bigratio) + (OldSolSvnC(i)*smallratio)
            AvgSolSvnLJ(i) = (AvgSolSvnLJ(i)*bigratio) + (OldSolSvnLJ(i)*smallratio)
            AvgSolSvnC2(i) = (AvgSolSvnC2(i)*bigratio) + (OldSolSvnC(i)**2*smallratio)
            AvgSolSvnLJ2(i) = (AvgSolSvnLJ2(i)*bigratio) + (OldSolSvnLJ(i)**2*smallratio)
            AvgSolSvnCF(i) = (AvgSolSvnCF(i)*bigratio) + (OldSolSvnCF(i)*smallratio)
            AvgSolSvnLJF(i) = (AvgSolSvnLJF(i)*bigratio) + (OldSolSvnLJF(i)*smallratio)
            AvgSolSvnCF2(i) = (AvgSolSvnCF2(i)*bigratio) + (OldSolSvnCF(i)**2*smallratio)
            AvgSolSvnLJF2(i) = (AvgSolSvnLJF2(i)*bigratio) + (OldSolSvnLJF(i)**2*smallratio)
            AvgSolSvnCB(i) = (AvgSolSvnCB(i)*bigratio) + (OldSolSvnCB(i)*smallratio)
            AvgSolSvnLJB(i) = (AvgSolSvnLJB(i)*bigratio) + (OldSolSvnLJB(i)*smallratio)
            AvgSolSvnCB2(i) = (AvgSolSvnCB2(i)*bigratio) + (OldSolSvnCB(i)**2*smallratio)
            AvgSolSvnLJB2(i) = (AvgSolSvnLJB2(i)*bigratio) + (OldSolSvnLJB(i)**2*smallratio)
          endif

          if (NGCSolutes.gt.0) then
            AvgSolGCSC(i) = (AvgSolGCSC(i)*bigratio) + (OldSolGCSC(i)*smallratio)
            AvgSolGCSLJ(i) = (AvgSolGCSLJ(i)*bigratio) + (OldSolGCSLJ(i)*smallratio)
            AvgSolGCSC2(i) = (AvgSolGCSC2(i)*bigratio) + (OldSolGCSC(i)**2*smallratio)
            AvgSolGCSLJ2(i) = (AvgSolGCSLJ2(i)*bigratio) + (OldSolGCSLJ(i)**2*smallratio)
            AvgSolGCSCF(i) = (AvgSolGCSCF(i)*bigratio) + (OldSolGCSCF(i)*smallratio)
            AvgSolGCSLJF(i) = (AvgSolGCSLJF(i)*bigratio) + (OldSolGCSLJF(i)*smallratio)
            AvgSolGCSCF2(i) = (AvgSolGCSCF2(i)*bigratio) + (OldSolGCSCF(i)**2*smallratio)
            AvgSolGCSLJF2(i) = (AvgSolGCSLJF2(i)*bigratio) + (OldSolGCSLJF(i)**2*smallratio)
            AvgSolGCSCB(i) = (AvgSolGCSCB(i)*bigratio) + (OldSolGCSCB(i)*smallratio)
            AvgSolGCSLJB(i) = (AvgSolGCSLJB(i)*bigratio) + (OldSolGCSLJB(i)*smallratio)
            AvgSolGCSCB2(i) = (AvgSolGCSCB2(i)*bigratio) + (OldSolGCSCB(i)**2*smallratio)
            AvgSolGCSLJB2(i) = (AvgSolGCSLJB2(i)*bigratio) + (OldSolGCSLJB(i)**2*smallratio)
         endif
        enddo
      endif                                    
      
c     now the solvent related averages
      if (NSolvents.gt.0) then
        AvgSvnSvnC = (AvgSvnSvnC*bigratio) + (OldSvnSvnC*smallratio)
        AvgSvnSvnLJ = (AvgSvnSvnLJ*bigratio) + (OldSvnSvnLJ*smallratio)
        AvgSvnSvnC2 = (AvgSvnSvnC2*bigratio) + (OldSvnSvnC**2*smallratio)
        AvgSvnSvnLJ2 = (AvgSvnSvnLJ2*bigratio) + (OldSvnSvnLJ**2*smallratio)
        AvgSvnSvnCF = (AvgSvnSvnCF*bigratio) + (OldSvnSvnCF*smallratio)
        AvgSvnSvnLJF = (AvgSvnSvnLJF*bigratio) + (OldSvnSvnLJF*smallratio)
        AvgSvnSvnCF2 = (AvgSvnSvnCF2*bigratio) + (OldSvnSvnCF**2*smallratio)
        AvgSvnSvnLJF2 = (AvgSvnSvnLJF2*bigratio) + (OldSvnSvnLJF**2*smallratio)
        AvgSvnSvnCB = (AvgSvnSvnCB*bigratio) + (OldSvnSvnCB*smallratio)
        AvgSvnSvnLJB = (AvgSvnSvnLJB*bigratio) + (OldSvnSvnLJB*smallratio)
        AvgSvnSvnCB2 = (AvgSvnSvnCB2*bigratio) + (OldSvnSvnCB**2*smallratio)
        AvgSvnSvnLJB2 = (AvgSvnSvnLJB2*bigratio) + (OldSvnSvnLJB**2*smallratio)
       
       if (NGCSolutes.gt.0) then
         do i=1,NSolvents
        AvgGCSSvnC(i) = (AvgGCSSvnC(i)*bigratio) + (OldGCSSvnC(i)*smallratio)
        AvgGCSSvnLJ(i) = (AvgGCSSvnLJ(i)*bigratio) + (OldGCSSvnLJ(i)*smallratio)
        AvgGCSSvnC2(i) = (AvgGCSSvnC2(i)*bigratio) + (OldGCSSvnC(i)**2*smallratio)
        AvgGCSSvnLJ2(i) = (AvgGCSSvnLJ2(i)*bigratio) + (OldGCSSvnLJ(i)**2*smallratio)
        AvgGCSSvnCF(i) = (AvgGCSSvnCF(i)*bigratio) + (OldGCSSvnCF(i)*smallratio)
        AvgGCSSvnCF2(i) = (AvgGCSSvnCF2(i)*bigratio) + (OldGCSSvnCF(i)**2*smallratio)
        AvgGCSSvnLJF(i) = (AvgGCSSvnLJF(i)*bigratio) + (OldGCSSvnLJF(i)*smallratio)
        AvgGCSSvnLJB(i) = (AvgGCSSvnLJB(i)*bigratio) + (OldGCSSvnLJB(i)*smallratio)
        AvgGCSSvnLJF2(i) = (AvgGCSSvnLJF2(i)*bigratio) + (OldGCSSvnLJF(i)**2*smallratio)
        AvgGCSSvnCB(i) = (AvgGCSSvnCB(i)*bigratio) + (OldGCSSvnCB(i)*smallratio)
        AvgGCSSvnCB2(i) = (AvgGCSSvnCB2(i)*bigratio) + (OldGCSSvnCB(i)**2*smallratio)
        AvgGCSSvnLJB2(i) = (AvgGCSSvnLJB2(i)*bigratio) + (OldGCSSvnLJB(i)**2*smallratio)
        
       
        enddo
       endif
      endif
      
c     now the gc-solute energy (gcs-gcs)
      if (NGCSolutes.gt.0) then
        AvgGCSGCSC = (AvgGCSGCSC*bigratio) + (OldGCSGCSC*smallratio)
        AvgGCSGCSLJ = (AvgGCSGCSLJ*bigratio) + (OldGCSGCSLJ*smallratio)
        AvgGCSGCSC2 = (AvgGCSGCSC2*bigratio) + (OldGCSGCSC**2*smallratio)
        AvgGCSGCSLJ2 = (AvgGCSGCSLJ2*bigratio) + (OldGCSGCSLJ**2*smallratio)
        AvgGCSGCSCF = (AvgGCSGCSCF*bigratio) + (OldGCSGCSCF*smallratio)
        AvgGCSGCSLJF = (AvgGCSGCSLJF*bigratio) + (OldGCSGCSLJF*smallratio)
        AvgGCSGCSCF2 = (AvgGCSGCSCF2*bigratio) + (OldGCSGCSCF**2*smallratio)
        AvgGCSGCSLJF2 = (AvgGCSGCSLJF2*bigratio) + (OldGCSGCSLJF**2*smallratio)
        AvgGCSGCSCB = (AvgGCSGCSCB*bigratio) + (OldGCSGCSCB*smallratio)
        AvgGCSGCSLJB = (AvgGCSGCSLJB*bigratio) + (OldGCSGCSLJB*smallratio)
        AvgGCSGCSCB2 = (AvgGCSGCSCB2*bigratio) + (OldGCSGCSCB**2*smallratio)
        AvgGCSGCSLJB2 = (AvgGCSGCSLJB2*bigratio) + (OldGCSGCSLJB**2*smallratio)
c
      endif
c     now call a routine to calculate the averages of the extra energies
      call calcAverageExtraEnergies(bigratio,smallratio)
                
      return
      end
