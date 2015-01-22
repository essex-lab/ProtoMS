      subroutine LREEnergy(lamB,lamF,nrg,nrgb,nrgf)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'enums.inc'
      include 'gb.inc'
c###########################################################
c
c     This subroutine calculates the energies needed
c     for lambda-replica exchange
c
c###########################################################

      double precision lamB,lamF
      double precision oldlamf,oldlamb
      double precision nrg,nrgb,nrgf
      double precision inrgb,inrgf
      double precision pme,pmeb,pmef
      double precision totalEnergy
      double precision ExtraEnergy

c     save the old lambdaf and lambdab
      oldlamf = lambdaf
      oldlamb = lambdab
            
c     set the new values of lambdaf and lambdab
      lambdaf = lamF
      lambdab = lamB
      
c     make sure that the values of lambda are sane!
      if (lambdaf.gt.ONE) lambdaf = ONE
      if (lambdaf.lt.ZERO) lambdaf = ZERO
      
      if (lambdab.gt.ONE) lambdab = ONE
      if (lambdab.lt.ZERO) lambdab = ZERO
      
c     now update the parameters and conformations of the system
      call recalculateParameters
      call recalculateGeometries
      
      call setAllMoved
c     The GB energy needs to be recalculated properly since it 
c     may have drifted
      if (doGB) then 
        call ZeroAllPDACont
        call ZeroPDACont
        call gbsavesystem
      endif
c     now calculate the total energy
c     ExtraEnergy must be set to NEW for the GB energy calculation 
c     to be correct. Note that it means that retienergy should be called
c     once at the end of a simulation, otherwise I am not too sure 
c     what this will cause to the running averages.
      nrg = totalEnergy(DONTCARE,nrgb,nrgf) + extraEnergy(NEW,inrgb,inrgf)
      nrgb = nrgb + inrgb
      nrgf = nrgf + inrgf

      call clearMoveList

c     PME
      if (pmeblocks.gt.0) then
        call pmePerturbeEnergy(pme,pmeb,pmef)
        nrg = nrg + pme
        nrgb = nrgb + pmeb
        nrgf = nrgf + pmef
      endif

c     now reset the parameters
      lambdaf = oldlamf
      lambdab = oldlamb
      call recalculateParameters
      call recalculateGeometries

      end
