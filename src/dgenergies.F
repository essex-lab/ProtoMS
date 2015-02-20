      subroutine dGEnergies
      include 'dimensions.inc'
      include 'printer.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'gb.inc'
      include 'averages.inc'
      include 'energies.inc'
      include 'metropolis.inc'
      include 'parallel.inc'
c###########################################################
c
c     This subroutine calculates energies for 
c     performing TI or BAR free energy estimates
c     It calculates the total energy for three different
c     lambda states (current,backward, and forward)
c     and calculate dV/dlambda using finite difference
c     prints it to the energy stream
c
c     (C) Samuel Genheden  20 February 2013
c
c###########################################################

      double precision oldlamf,oldlamb,oldlam
      double precision lb,lf,dvdl,dgf,dgb,dvdl2,dl
      double precision nrg,nrgb,nrgf
      double precision inrgb,inrgf
      double precision nrgti,nrgtib,nrgtif
      double precision totalEnergy
      double precision pme,pmeb,pmef
      double precision ExtraEnergy,LambdaDerivative
      integer lambindex1, lambindex2,i

      if (.not.(doLambdaRE.and.printFE)) return

      if (mbar) then
        lambindex1 = 1
        lambindex2 = nlambdas
      else
c       setting lambda backwards index
        if (mylambda.gt.1) then 
          lambindex1 = mylambda-1
        else
          lambindex1 = mylambda
        endif
c       setting lambda forwards index
        if (mylambda.lt.nlambdas) then
          lambindex2 = mylambda+1
        else
          lambindex2 = mylambda
        endif
      endif

      if (reCalculateFreeEnergies) then
        reCalculateFreeEnergies = .false.
        
c       save the old lambdaf and lambdab
        oldlamf = lambdaf
        oldlam = lambda
        oldlamb = lambdab

c       doing the gradient calculations
        dl = max(lambdaf-lambda,lambda-lambdab)
        dgf=-dlog(AVGFORWARDSENERGY) / mt_beta_T
        dgb=-dlog(AVGBACKENERGY) / mt_beta_T
        if (lambda.eq.ZERO) then
          numericalGradient = dgf / dl
        else if (lambda.eq.ONE) then
          numericalGradient = -dgb / dl
        else
          numericalGradient = (dgf-dgb)/(2.0d0*dl)
        endif   
        AnalyticalGradient = lambdaderivative() 


c       start the bar / mbar calculations
        do i=lambindex1,lambindex2,3

          lambdaf = lambdaval(min(i+2,lambindex2))
          lambda = lambdaval(min(i+1,lambindex2))
          lambdab = lambdaval(i)

          if (lambdaf.gt.ONE) lambdaf = ONE
          if (lambdaf.lt.ZERO) lambdaf = ZERO
      
          if (lambdab.gt.ONE) lambdab = ONE
          if (lambdab.lt.ZERO) lambdab = ZERO
      
          call recalculateParameters
          call recalculateGeometries

          call setAllMoved
          if (doGB) then 
            call ZeroAllPDACont
            call ZeroPDACont
            call gbsavesystem
          endif
          nrg = totalEnergy(DONTCARE,nrgb,nrgf) + extraEnergy(NEW,inrgb,inrgf)
          nrgb = nrgb + inrgb
          nrgf = nrgf + inrgf

c         PME
          if (pmeblocks.gt.0) then
            call pmePerturbeEnergy(pme,pmeb,pmef)
            nrg = nrg + pme
            nrgb = nrgb + pmeb
            nrgf = nrgf + pmef
          endif

          MBAREnergy(min(i+2,lambindex2)) = nrgf
          MBAREnergy(min(i+1,lambindex2)) = nrg
          MBAREnergy(i) = nrgb
        enddo
      
        call clearMoveList
c       now reset the parameters
        lambdaf = oldlamf
        lambda = oldlam
        lambdab = oldlamb
        call recalculateParameters
        call recalculateGeometries

c     Now we have up to date energies to print
      endif   

c     write headers
      call prettyLine(RESULTS,'FREE ENERGY DATA')
      call prettyLine(RESULTS,'lambda | Total Energy')
      
      do i=lambindex1,lambindex2
        write(printstring,'(F5.3,X,A,X,F12.4)') lambdaval(i),'|',MBAREnergy(i)
        call printLine(RESULTS,printstring)
      enddo
      
      write(printstring,'(A,X,F12.4,X,F12.4)') 'Gradient',NumericalGradient,AnalyticalGradient
      call prettyline(RESULTS,printstring)

      end
