      double precision function calcJPMF(isol,lambda)
      include 'dimensions.inc'
      include 'solutes.inc'

c###############################################################
c
c     This function calculate the PMF solvation bias
c     applied to a theta move of a JAWS-1 simulation of solutes
c     to account for solvation effects. lambda = 1 - theta
c
c     (C) Samuel Genheden  October 2012
c
c###############################################################      
      

      integer isol,i
      double precision lambda,lambdamul,lambda_real

      if (SolJPmfOrd(isol).eq.0) then
        calcJPMF = 0.0d0
        return
      endif
     
      calcJPMF = SolJPmfCoeff(isol,1)
      lambdamul = lambda
      do i=2,SolJPmfOrd(isol)
        calcJPMF = calcJPMF + lambdamul*SolJPmfCoeff(isol,i)
        lambdamul = lambdamul*lambda
      enddo
      return
     
      end
