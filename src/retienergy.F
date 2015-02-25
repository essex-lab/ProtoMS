      subroutine retiEnergy(retideltau,retideltal)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'enums.inc'
      include 'gb.inc'
c###########################################################
c
c     This subroutine calculates the energies needed
c     for RETI. It uses the value of 'retidelta' to
c     change lambdaf and lambdab to lambda+retidelta
c     and lambda-retidelta respectively. It then updates
c     the conformations to the new lambda values, and 
c     calculates the new total energy. It prints the result
c     to the RETI stream
c
c     It then resets the values of lambda back to where
c     they were
c
c     (C) Christopher Woods  5 February 2002
c     Modified by J Michel November 2006
c
c     The retidelta increment is not necessary the same for the forward 
c     and backward simulations
c     the new syntax is call retienergy retideltau retidelal
c###########################################################

      double precision oldlamf,oldlamb
      double precision retideltau,retideltal
      double precision nrg,nrgb,nrgf,pme,pmeb,pmef
      double precision inrgb,inrgf
      double precision totalEnergy
      double precision ExtraEnergy

c     save the old lambdaf and lambdab
      oldlamf = lambdaf
      oldlamb = lambdab
      
      write(printstring,*) "Calculating RETI energies with delta = ",retideltau,retideltal
      call printLine(INFO,printstring)
      
c     set the new values of lambdaf and lambdab
      lambdaf = lambda + retideltau
      lambdab = lambda - retideltal
      
c     make sure that the values of lambda are sane!
      if (lambdaf.gt.ONE) lambdaf = ONE
      if (lambdaf.lt.ZERO) lambdaf = ZERO
      
      if (lambdab.gt.ONE) lambdab = ONE
      if (lambdab.lt.ZERO) lambdab = ZERO
      
c     now update the parameters and conformations of the system
      call recalculateParameters
      call recalculateGeometries

      
      call setAllMoved
c     The GB energy needs to be recalculated properly since it 
c     may have drifted
      if (doGB) then 
        call ZeroAllPDACont
        call ZeroPDACont
        call gbsavesystem
      endif
c     now calculate the total energy
c     ExtraEnergy must be set to NEW for the GB energy calculation 
c     to be correct. Note that it means that retienergy should be called
c     once at the end of a simulation, otherwise I am not too sure 
c     what this will cause to the running averages.
      nrg = totalEnergy(DONTCARE,nrgb,nrgf) + extraEnergy(NEW,inrgb,inrgf)
      nrgb = nrgb + inrgb
      nrgf = nrgf + inrgf

      call clearMoveList

c     PME
      if (pmeblocks.gt.0) then
        call pmePerturbeEnergy(pme,pmeb,pmef)
        nrg = nrg + pme
        nrgb = nrgb + pmeb
        nrgf = nrgf + pmef
      endif

c     now write this all out to the RETI stream
      write(printstring,'(A,6(F20.8,X))') "RETINRG ",lambda,nrg,
     .        lambdab,nrgb,lambdaf,nrgf
      call printLine(RETI,printstring)
      call printLine(INFO,printstring)
c      call printLine(ENERGY,printstring)
c     now reset the parameters
      lambdaf = oldlamf
      lambdab = oldlamb
      call recalculateParameters
      call recalculateGeometries

      end
