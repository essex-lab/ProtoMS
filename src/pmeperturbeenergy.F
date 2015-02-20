      subroutine pmePerturbeEnergy(nrg,nrgb,nrgf)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'enums.inc'
      include 'averages.inc'
c###########################################################
c
c     This routine performs PME calculations
c     for a perturbation, i.e. for backward
c     and forward as well as the current lambda
c     This should only be done once, as it is very costly
c
c     Samuel Genheden, November 2013
c
c###########################################################

      double precision oldlamf,oldlamb,oldlam
      double precision retideltau,retideltal
      double precision nrg,nrgb,nrgf
      double precision selfnrg,exclnrg,meshnrg,selfsvnnrg,exclsvnnrg,chnrg,pmeEnergy

c     c
c     Calculate for the current lambda
c     c
      nrg = pmeEnergy(selfnrg,exclnrg,meshnrg,chnrg,selfsvnnrg,exclsvnnrg,.true.)

      write(printstring,11)'PME-energy',nrg,selfnrg,exclnrg,meshnrg,chnrg
 11   format(A20,5(F15.4,X))
      call printLine(ENERGY,printstring)
      call printLine(INFO,printstring)


c     save the old lambdaf and lambdab and lambda
      oldlamf = lambdaf
      oldlamb = lambdab
      oldlam  = lambda

c     c
c     Start doing the backwards lambda            
c     c

c     set the new values of lambdaf and lambdab and lambda
      lambdaf = lambdab
      lambda  = lambdab
           
c     now update the parameters and conformations of the system
      call recalculateParameters
      call recalculateGeometries    
      call setAllMoved

      nrgb = pmeEnergy(selfnrg,exclnrg,meshnrg,chnrg,selfsvnnrg,exclsvnnrg,.true.)

      write(printstring,11)'PME-backwards',nrgb,selfnrg,exclnrg,meshnrg,chnrg
      call printLine(ENERGY,printstring)
      call printLine(INFO,printstring)


c     c
c     Start doing the forwards lambda
c     c

      lambdaf = oldlamf
      lambdab = oldlamf
      lambda  = oldlamf

c     now update the parameters and conformations of the system
      call recalculateParameters
      call recalculateGeometries    
      call setAllMoved

      nrgf = pmeEnergy(selfnrg,exclnrg,meshnrg,chnrg,selfsvnnrg,exclsvnnrg,.true.)

      write(printstring,11)'PME-forwards',nrgf,selfnrg,exclnrg,meshnrg,chnrg
      call printLine(ENERGY,printstring)
      call printLine(INFO,printstring)


c     now reset the parameters
      lambdaf = oldlamf
      lambdab = oldlamb
      lambda  = oldlam

      call recalculateParameters
      call recalculateGeometries

      write(printstring,14)'PME-total',nrg,nrgb,nrgf
 14   format(A20,3(F15.4,X))
      call printLine(ENERGY,printstring)
      call printLine(INFO,printstring)

      AVGPMEENERGY = nrg
      AVGPMEENERGYb = nrgb
      AVGPMEENERGYf = nrgf

      end
