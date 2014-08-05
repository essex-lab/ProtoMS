      subroutine dGEnergies(mbar)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'enums.inc'
      include 'gb.inc'
      include 'averages.inc'
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
      double precision ExtraEnergy,LambdaDerivative
      integer lambindex1, lambindex2,i
      logical mbar


c     save the old lambdaf and lambdab
      oldlamf = lambdaf
      oldlam = lambda
      oldlamb = lambdab

c     doing the gradient calculations
c     we will print them at the very end
      nrgti = totalEnergy(DONTCARE,nrgtib,nrgtif) + extraEnergy(DONTCARE,inrgb,inrgf)
      nrgtib = nrgtib + inrgb
      nrgtif = nrgtif + inrgf
      dl = max(lambdaf-lambda,lambda-lambdab)
      dgf=-dlog(AVGFORWARDSENERGY) / mt_beta_T
      dgb=-dlog(AVGBACKENERGY) / mt_beta_T
      if (lambda.eq.ZERO) then
        dvdl = (nrgtif - nrgti) / dl
        dvdl2 = dgf / dl
      else if (lambda.eq.ONE) then
        dvdl = (nrgti - nrgtib) / dl
        dvdl2 = -dgb / dl
      else
        dvdl = (nrgtif - nrgtib) / (2.0*dl)
        dvdl2 = (dgf-dgb)/(2.0d0*dl)
      endif    


c     write headers
      call prettyLine(RESULTS,'FREE ENERGY DATA')
      call prettyLine(RESULTS,'lambda | Total Energy')


c     start the bar / mbar calculations
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
c     do the actual perturbations
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

        write(printstring,'(F5.3,X,A,X,F12.4)') lambdab,'|',nrgb
        call printLi ne(RESULTS,printstring)
        if (.not.(lambda.eq.lambdab)) then
          write(printstring,'(F5.3,X,A,X,F12.4)') lambda,'|',nrg
          call printLine(RESULTS,printstring)
        endif
        if (.not.(lambdaf.eq.lambda)) then
          write(printstring,'(F5.3,X,A,X,F12.4)') lambdaf,'|',nrgf
          call printLine(RESULTS,printstring)
        endif

      enddo
      
      call clearMoveList
c     now reset the parameters
      lambdaf = oldlamf
      lambda = oldlam
      lambdab = oldlamb
      call recalculateParameters
      call recalculateGeometries

c      call prettyLine(RESULTS,'Gradient. Instantaneous | On running averages')
c      write(printstring,'(F12.4,X,A,X,F12.4)') dvdl,'|',dvdl2
      write(printstring,'(A,X,F12.4,X,F12.4)') 'Gradient',dvdl2,lambdaderivative()
      call prettyline(RESULTS,printstring)

      end
