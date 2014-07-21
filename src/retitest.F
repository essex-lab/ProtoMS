      logical function retiTest(lamA,lamfA,lambA,nrgA,nrgfA,nrgbA,lamB,lamfB,lambB,nrgB,nrgfB,nrgbB)
      include 'dimensions.inc'
      include 'metropolis.inc'
      include 'constants.inc'
      include 'printer.inc'
c###################################################################
c
c     This routine performs a Metropolis Test to see if two
c     Hamiltonians should be swapped
c
c     Samuel Genheden, 2012
c
c###################################################################    
      double precision lamA,lamfA,lambA,nrgA,nrgfA,nrgbA,lamB,lamfB,lambB,nrgB,nrgfB,nrgbB
      double precision nrgAi,nrgAj,nrgBi,nrgBj,expval,ranval,delta
      double precision thres,grnd
      parameter (thres=0.00001)
      
      retiTest = .false.
      
      nrgAi = nrgA
      if (abs(lamA-lambB).lt.thres) then
        nrgAj = nrgbB
      else if (abs(lamA-lamfB).lt.thres) then
        nrgAj = nrgfB
      else
        write(printstring,"(A)")"Could not determine nrgAj. Something is seriously wrong! Exit."
        call prettyLine(FATAL,printstring)     
        return
      endif
      
      nrgBj = nrgB
      if (abs(lamB-lambA).lt.thres) then
        nrgBi = nrgbA
      else if (abs(lamB-lamfA).lt.thres) then
        nrgBi = nrgfA
      else
        write(printstring,"(A)")"Could not determine nrgBi. Something is seriously wrong! Exit."
        call prettyLine(FATAL,printstring)     
        return
      endif
      
      delta = -((nrgBj-nrgBi)+(nrgAi-nrgAj))
      if (delta.lt.ZERO)  then
        retiTest = .true.
        return
      endif

      ranval = grnd()
      expval = dexp(-mt_beta_T*delta)
c      write(printstring,'(3F4.2,3F12.3,3F4.2,3F12.3,2F10.3)')lamA,lamfA,lambA,nrgA,nrgfA,nrgbA,lamB,lamfB,lambB,nrgB,nrgfB,nrgbB,expval,ranval
c      call prettyline(INFO,printstring)

c      write(6,*)"E",ranval,expval,delta,mt_beta_T     
      if (expval.ge.ranval) then
        retiTest = .true.
      endif
      
      return
      
      end function

