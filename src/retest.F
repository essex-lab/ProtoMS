      logical function reTest(betaA,betaB,eneA,eneB)
      include 'dimensions.inc'
      include 'metropolis.inc'
c      include 'reti.inc'
      include 'constants.inc'
      include 'printer.inc'
c###################################################################
c
c     This routine performs a Metropolis Test to see if two
c     Hamiltonians should be swapped
c
c     Samuel Genheden, 2013
c
c###################################################################    
      double precision betaA,betaB,eneA,eneB
      double precision expval,ranval,grnd,test

      test = (betaB-betaA)*(eneA-eneB)   

c      write(6,*)" eneA ",eneA," eneB ",eneB


      reTest = .false.
      ranval = grnd()
      expval = dexp(-test)
      if (expval.ge.ranval) then
        reTest = .true.
      endif

      
      return

      end function

