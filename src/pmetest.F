      subroutine pmeTest(oldnrg,newnrg,selfsvnnrg,exclsvnnrg,block)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'simparams.inc'
      include 'solvents.inc'
      include 'moves.inc'
      include 'metropolis.inc'
      include 'enums.inc'
      include 'energies.inc'
c###########################################################
c
c     This function performs PME calculations and
c     Metropolis test to see if the new configuration should
c     be accepted
c
c     (C) Samuel Genheden        October 2013
c
c###########################################################
      integer block
      double precision oldnrg,newnrg,tmp
      double precision selfnrg,exclnrg,meshnrg,chnrg,selfsvnnrg,exclsvnnrg,pmeEnergy

      character ctmp
      logical templ,metropolis

      if (block.eq.0) then
c       This is what we will compare to in the first Metropolis test
        oldnrg = pmeEnergy(selfnrg,exclnrg,meshnrg,chnrg,selfsvnnrg,exclsvnnrg,.true.)
        write(printstring,11)0,OLDTOTALENERGY,selfnrg,exclnrg,meshnrg,chnrg,oldnrg,'A'
 11     format('PME ',I5,X,6(F15.4,X),X,A)
        call printLine(ENERGY,printstring)
        
c       In case the Metropolis test fails, we need to save the old configuration
        call saveState()

      else
c       Calculate the Ewald (PME) energy of the system
        newnrg = pmeEnergy(selfnrg,exclnrg,meshnrg,chnrg,selfsvnnrg,exclsvnnrg,.false.)        
c       Do the Metropolis test
        tmp = mt_bias
        mt_bias = 1.0d0
        templ = metropolis(oldnrg,newnrg,0)
        mt_bias = tmp
        ctmp = 'F'
        if (templ) ctmp = 'A'
        write(printstring,11)block,OLDTOTALENERGY,selfnrg,exclnrg,meshnrg,chnrg,newnrg,ctmp
        call printLine(ENERGY,printstring)
        NEwaldTry = NEwaldTry + 1
          
c       If we succeed we are just carrying on
        if (templ) then
          oldnrg = newnrg
          NEwaldAcp = NEwaldAcp + 1
c         Write current configuration
          call saveState()
c       But if we fail, we need to restore the old configuration and re-calculate total energy
        else
          call prettyLine(INFO,"Failed Ewald test")
          call loadState()
        endif
      
      endif

      return

      end subroutine
