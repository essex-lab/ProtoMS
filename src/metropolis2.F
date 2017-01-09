      logical function metropolis2(olde,ienew,movtyp)
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
      include 'flags.inc'
      include 'moves.inc'
      include 'energies.inc'
      include 'stacks.inc'

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
c     movtyp == 2        Addition for GC test 
c     movtyp == 3        Deletion for GC test      
c
c##########################################################

      double precision ienew,delta,expval
      double precision ranval,grnd,olde  
      integer movtyp

      metropolis2 = .false.

       if (movtyp.le.0) then
c#########################################################
c######## NVT metropolis MC acceptance test ##############
        delta = ienew - olde
        

c       get the random number
        ranval = grnd()

        expval = ((gridx*gridy*gridz/vavail)*(dexp(B))/GCSCavity)*dexp((olde-ienew)*mt_beta_T)

        if (expval.lt.ranval) then
c         move failed
          metropolis2 = .false.
          mt_bias = ONE
          return
        else
c         move passed
          metropolis2 = .true.
          mt_bias = ONE
          return
        endif
        
      endif
      
      mt_bias = ONE
      
      return
      end
