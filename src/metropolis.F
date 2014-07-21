      logical function metropolis(eold,enew,movtyp)
      include 'dimensions.inc'
      include 'metropolis.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'solvents.inc'
      include 'enums.inc'
      include 'pressure.inc'
      include 'simparams.inc'
c##########################################################
c
c     This function implements the standard metropolis
c     test, and returns whether it passes or fails.
c
c     The old energy is in eold, the new energy is
c     in enew, and 'movtyp' is an integer giving the
c     type of move (a placeholder for now that may
c     allow different tests to be implemented, e.g.
c     NPT metropolis, multicanonical etc...
c
c     There is a common block (metropolis.inc) that 
c     contains additional space for variables that
c     may be used in any tests (e.g. pressures or volumes)
c
c     In this way, variables for tests can be added without changing
c     the method of calling this function.
c
c     movtyp <= 0        Normal NVT metropolis test (load up mt_beta_T)
c     movtyp == 1        Normal NPT met test (load mt_vold and mtv_new
c     movtyp == 2        Multicanonical test (not yet implemented...)
c
c##########################################################

      double precision eold,enew,delta,expval
      double precision ranval,grnd,PMFold,PMFnew
      double precision idelv,ivratio
      integer movtyp

      metropolis = .false.

      if (movtyp.le.0) then
c#########################################################
c######## NVT metropolis MC acceptance test ##############
c In this implementation, we now add on the PMF information. Here, it's
c for water...

       delta = enew - eold
       
       if (delta.le.ZERO) then
c         move is automatically accepted
          metropolis = .true.
          mt_bias = ONE
          return
        endif
c             PRINT*, mt_beta_T
c       get the random number
        ranval = grnd()
c       compare exp(-dE / kT) to random number between 0 and 1
        expval = mt_bias*dexp(-mt_beta_T*delta)
c        write(printstring,*) 'expval ',expval,' ranval ',ranval, 'B ',B, -mt_beta_T
c        call printLine(DEBUG,printstring)
        if (expval.lt.ranval) then
c         move failed
          metropolis = .false.
          mt_bias = ONE
          return
        else
c         move passed
          metropolis = .true.
          mt_bias = ONE
          return
        endif
        
      else if (movtyp.eq.NPTMCTEST) then
c####################################################
c######### NPT metropolis MC test ###################
        delta = enew-eold
        idelv = intPressure*(mt_vnew - mt_vold)
        ivratio = dble(NSolvents+NProteins+NSolutes)*dlog(mt_vnew/mt_vold)
        
        expval = mt_bias*dexp( -mt_beta_T*(delta+idelv) + ivratio)
        
        if (expval.gt.ONE) then
c         move automatically passed
          metropolis = .true.
          mt_bias = ONE
          return
        endif
        
        ranval = grnd()
        
        if (expval.lt.ranval) then
c         move failed
          metropolis = .false.
          mt_bias = ONE
          return
        else
c         move passed
          metropolis = .true.
          mt_bias = ONE
          return
        endif
        
      else
c####################################################
c###### Requested test not recognised ###############
        write(printstring,*) "This type of MC test has not yet been implemented...",movtyp
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is probably a programming bug!")
        metropolis = .false.
      endif
      
      mt_bias = ONE
      
      return
      
      end
