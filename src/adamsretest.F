      logical function adamsreTest(adamsA,adamsB,ngcsolsA,ngcsolsB)
      include 'dimensions.inc'
      include 'metropolis.inc'
c      include 'reti.inc'
      include 'constants.inc'
      include 'printer.inc'
c###################################################################
c
c     This routine performs a Metropolis Test to see if two
c     GC simulations with Adams potentials should be swapped.
c
c     Gregory Ross, 2015
c	  Modified from retest.F	
c###################################################################    
      double precision ngcsolsA,ngcsolsB,adamsA,adamsB,beta
      double precision expval,ranval,grnd,test

      test = (ngcsolsA-ngcsolsB)*(adamsB-adamsA)*mt_beta_T   

      adamsretest = .false.
      ranval = grnd()
      expval = dexp(test)
      if (expval.ge.ranval) then
        adamsretest = .true.
      endif

      
      return

      end function

