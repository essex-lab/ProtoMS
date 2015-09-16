      double precision function calcJPMF(isol,theta)
      include 'dimensions.inc'
      include 'solutes.inc'
      include 'flags.inc'
      include 'simparams.inc'

c###############################################################
c
c     This function calculate the PMF solvation bias
c     applied to a theta move of a JAWS-1 simulation of solutes
c     to account for solvation effects
c
c     (C) Samuel Genheden  October 2012
c
c###############################################################      
      

      integer isol,i
      double precision theta,thetamul,theta_real

      if (SolJPmfOrd(isol).eq.0) then
        calcJPMF = 0.0d0
        return
      endif
     
      calcJPMF = SolJPmfCoeff(isol,1)
      thetamul = theta
      do i=2,SolJPmfOrd(isol)
        calcJPMF = calcJPMF + thetamul*SolJPmfCoeff(isol,i)
        thetamul = thetamul*theta
      enddo
      if (doREST) then
        calcJMPF = calcJPMF * restscaling(7)
      endif
      return
     
      end
